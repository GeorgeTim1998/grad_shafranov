# В этой задаче решаем динамическую задачу при наличии точечных управляющих катушек, камеры с магнитной проницаемостью отличной от нуля


# %% Imports
from fenics import *
import matplotlib.pyplot as matplt
import logger
import mshr
import time
import funcs as fu
import MEPHIST_data as M
import logger
from geometry import Geometry
from boundary_conditions import BoundaryConditions
import point_source_data as psd
import numpy
from expressions import Expressions
from ITER_params import Problem

# %% Pre-programm stuff
t0 = time.time()
current_pyfile = "\n\n---------ITER.py---------"
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'ITER'

# %% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()
p = Problem()
e = Expressions()
fu.What_time_is_it(t0, 'Problem params logged')

# %% Domain and mesh definition
domain = geometry.rectangle_domain(
    area=[p.domain_geometry[0],
          p.domain_geometry[1],
          p.domain_geometry[2],
          p.domain_geometry[3]])

geometry.register_plot_domain(p.plot_domain)
iter_inner_surface = geometry.iter_FW()
iter_outer_surface = geometry.outer_iter_vessel()

vacuum_vessel = iter_outer_surface - iter_inner_surface

domain.set_subdomain(1, vacuum_vessel)  # vessel
domain.set_subdomain(2, iter_inner_surface)  # plasma

geometry.generate_mesh_in_domain(domain=domain, density=p.mesh_density)

# fu.fenics_plot(p, geometry.mesh, PATH, limits=1)

markers = MeshFunction("size_t", geometry.mesh,
                       geometry.mesh.topology().dim(), geometry.mesh.domains())

fu.fenics_plot(p, markers, PATH)
fu.fenics_plot(p, markers, "%s_nobar" % PATH)
fu.What_time_is_it(t0, 'Markers of domains plotted')

# %% Step coefficients classes


class Permeability(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers

    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 0:
            values[0] = p.VACUUM_PERMEABILITY  # vacuum
        elif self.markers[cell.index] == 1:
            values[0] = p.VESSEL_PERMEABILITY  # vessel
        else:
            values[0] = p.PLASMA_PERMEABILITY  # plasma

    def value_shape(self):
        return ()


class Conductivity(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers

    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 0:
            values[0] = p.VACUUM_CONDUCTIVITY  # vacuum
        elif self.markers[cell.index] == 1:
            values[0] = p.VESSEL_CONDUCTIVITY  # vessel
        else:
            values[0] = p.PLASMA_CONDUCTIVITY  # plasma

    def value_shape(self):
        return ()


# %% Define function space and step coefficients
V = FunctionSpace(geometry.mesh, 'P', 1)  # standard triangular mesh

mu = Permeability(geometry.mesh, degree=0)
sg = Conductivity(geometry.mesh, degree=0)

u = TrialFunction(V)  # u must be defined as function before expression def
v = TestFunction(V)

# %% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary)

# %% Solve
[r_2, r] = geometry.operator_weights(V)

dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

# %% Solve stationary
source = e.iter_point_source(problem=p, a=[p.disp_x[0], p.disp_z[0]])
a = dot(grad(u)/r, grad(r_2*v))*dx
L = source*r*v*dx(2)  # !!!

u0 = Function(V)
solve(a == L, u0, bc)
p.find_levels(u0, step=p.step,
              shrink=p.shrink,
              shrink_step=p.shrink_step)

fu.What_time_is_it(t0, 'Initial problem solved')
fu.countour_plot_via_mesh(geometry, u0, levels=p.levels,
                          PATH=PATH,
                          current_disp=p.centre_point,
                          plt_vessel=True,
                          do_plasma_centre=True,
                          colorbar=True,
                          xticks_array=p.xticks,
                          grid=True)
fu.What_time_is_it(t0, 'Plot with bar plotted')
fu.countour_plot_via_mesh(geometry, u0, levels=p.levels,
                          PATH=PATH+'_nobar',
                          current_disp=p.centre_point,
                          plt_vessel=True,
                          do_plasma_centre=True,
                          colorbar=False,
                          xticks_array=p.xticks,
                          grid=True)
fu.What_time_is_it(t0, 'Plot with no bar plotted')

dt = numpy.diff(p.t)
for i in range(len(dt)):
    logger.info("i = %d out of %d" % (i, len(dt)))
    fu.print_colored_n_white(colored_text="Time: ",
                             color='blue', white_text=str(p.t[i+1]))
    u = Function(V)
    v = TestFunction(V)

    source = e.iter_point_source(problem=p,
                                 a=[p.disp_x[i+1], p.disp_z[i+1]])

    # F = dot(grad(u)/r, grad(r_2*v))*dx \
    #     + fu.M0*mu*sg / dt[i] * (u - u0)*r*v*dx \
    #     - source*r*v*dx(2)
    F = dot(grad(u)/r, grad(r_2*v))*dx \
        - source*r*v*dx(2)

    solve(F == 0, u, bc)
    # p.find_levels(u, step=p.step)

    fu.What_time_is_it(t0, "Problem solved for t = %f" % p.t[i+1])
    # p.find_levels_with_exponent(u, exponent=p.exponent, step=p.step)
    fu.countour_plot_via_mesh(geometry, u, levels=p.levels,
                              PATH=PATH+'_nobar',
                              current_disp=[p.disp_x[i+1], p.disp_z[i+1]],
                              plt_vessel=True,
                              do_plasma_centre=True,
                              colorbar=False,
                              grid=True,
                              xticks_array=p.xticks)
    fu.countour_plot_via_mesh(geometry, u, levels=p.levels,
                              PATH=PATH,
                              current_disp=[p.disp_x[i+1], p.disp_z[i+1]],
                              plt_vessel=True,
                              do_plasma_centre=True,
                              colorbar=True,
                              grid=True,
                              xticks_array=p.xticks)

    u0 = u

    fu.print_colored(text="Iteration finished. %d/%d = %.2f" 
                     % (i, len(dt)-1, 100*i/(len(dt)-1)), color='yellow')


fu.What_time_is_it(t0, message='Done')
logger.log_n_output_colored_message(
    colored_message="'Done'\n", color='red', white_message='')

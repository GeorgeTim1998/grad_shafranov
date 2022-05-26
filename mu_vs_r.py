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
from mu_vs_r_problem_params import Problem

# %% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------mu_vs_r.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'mu_vs_r'

# %% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()
p = Problem()

# %% Domain and mesh definition
domain = geometry.rectangle_domain(
    area=[p.domain_geometry[0], p.domain_geometry[1], p.domain_geometry[2], p.domain_geometry[3]])
circle1 = geometry.circle_domain(centre_point=p.centre_point,
                                 radius=p.ves_inner_radius, segments=p.ves_inner_segments)
circle2 = geometry.circle_domain(centre_point=p.centre_point,
                                 radius=p.ves_outer_radius, segments=p.ves_outer_segments)

ring = circle2 - circle1

domain.set_subdomain(1, ring)  # vessel
domain.set_subdomain(2, circle1)  # plasma

geometry.generate_mesh_in_domain(domain=domain, density=p.mesh_density)

# fu.fenics_plot(geometry.mesh, PATH, '', '')

markers = MeshFunction("size_t", geometry.mesh,
                       geometry.mesh.topology().dim(), geometry.mesh.domains())

fu.fenics_plot(markers, PATH, '', '')
fu.fenics_plot(markers, "%s_nobar" % PATH, '', '')

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

# fu.fenics_plot(interpolate(sg, V), PATH, '', 'colorbar')
# fu.countour_plot_via_mesh(geometry, interpolate(
#     mu, V), levels=p.levels, PATH=PATH, plot_title='Permeability')
# fu.countour_plot_via_mesh(geometry, interpolate(
#     sg, V), levels=p.levels, PATH=PATH, plot_title='Conductivity')

u = TrialFunction(V)  # u must be defined as function before expression def
v = TestFunction(V)

# %% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary)

# %% Solve
[r_2, r] = geometry.operator_weights(V)

boundary_conditions.spheromak_boundary_condition(
    psi_0=p.psi_0, R=p.R, alpha=p.alpha)
point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(1)))

dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

# %% Solve stationary
a = dot(grad(u)/r, grad(r_2*v))*dx
L = sum(point_sources[2:len(point_sources)])*r*v*dx(0) + \
    boundary_conditions.spheromak_right_hand_expr*r*v*dx(2)

u0 = Function(V)
solve(a == L, u0, bc)

fu.What_time_is_it(t0, 'Initial problem solved')
fu.countour_plot_via_mesh(geometry, u0, levels=p.levels,
                          PATH=PATH, plot_title='')

dt = numpy.diff(p.t)
for i in range(len(dt)):
    u = Function(V)
    v = TestFunction(V)

    F = dot(grad(u)/r, grad(r_2*v))*dx + \
        fu.M0*mu*sg / dt[i] * (u - u0)*r*v*dx - \
        sum(point_sources[2:len(point_sources)])*r*v*dx(0) - \
        boundary_conditions.spheromak_right_hand_expr/p.R_t[i+1]*r*v*dx(2)

    solve(F == 0, u, bc)

    fu.What_time_is_it(t0, "Problem solved for t = %f" % p.t[i+1])
    fu.countour_plot_via_mesh(geometry, u, levels=p.levels,
                              PATH=PATH, plot_title='')
    fu.countour_plot_via_mesh_nocolorbar(geometry, u, levels=p.levels,
                                         PATH="%s_nobar" % PATH, plot_title='')

    u0 = u


logger.log_n_output_colored_message(
    colored_message="'Done'\n", color='red', white_message='')

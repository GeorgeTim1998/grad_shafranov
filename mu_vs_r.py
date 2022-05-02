#%% Imports
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
    
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------mu_vs_r.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'mu_vs_r'

#%% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()

VACUUM_PERMEABILITY = 1
VESSEL_PERMEABILITY = 5
logger.log_n_output_colored_message(colored_message="VACUUM_PERMEABILITY = ", color='green', white_message=str(VACUUM_PERMEABILITY))
logger.log_n_output_colored_message(colored_message="VESSEL_PERMEABILITY = ", color='green', white_message=str(VESSEL_PERMEABILITY))

levels = 20

#%% Domain and mesh definition
geometry.rectangle_mesh_init(r1 = 0.05, r2 = 0.55, z1 = -0.4, z2 = 0.4, default_mesh = 100) # Howell2014

domain = geometry.rectangle_domain(area=[0.05, 0.55, -0.4, 0.4])
circle1 = geometry.circle_domain(radius=0.1, segments=60)
circle2 = geometry.circle_domain(radius=0.15, segments=100)

ring = circle2 - circle1

domain.set_subdomain(1, ring)
domain.set_subdomain(2, circle1)

geometry.generate_mesh_in_domain(domain=domain, density=128)

plot(geometry.mesh)
fu.save_contour_plot(PATH, '')

markers = MeshFunction("size_t", geometry.mesh, geometry.mesh.topology().dim(), geometry.mesh.domains())

plot(markers)
fu.save_contour_plot(PATH, '')

class Permeability(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 1:
            values[0] = VESSEL_PERMEABILITY # vessel
        else:
            values[0] = VACUUM_PERMEABILITY # vacuum

class StepFunction(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 1:
            values[0] = VESSEL_PERMEABILITY # vessel
        else:
            values[0] = VACUUM_PERMEABILITY # vacuum

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh

mu = Permeability(geometry.mesh, degree=0)
fu.countour_plot_via_mesh(geometry, interpolate(mu, V), levels = levels, PATH = PATH, plot_title = '')

u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

# dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(1)))
[p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_pow_2(p_correction=100, F_correction=1)

L = mu*sum(point_sources)*r*v*dx 

a = dot(grad(u)/r, grad(r_2*v))*dx - mu*(p_coeff*r*r + F_2_coeff)*u*r*v*dx
u = Function(V)
solve(a == L, u, bc)

#%% Post solve
fu.What_time_is_it(t0, 'Variational problem solved')
fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')

fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
logger.info("'Done'"+"\n")
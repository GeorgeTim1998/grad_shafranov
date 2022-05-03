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
current_pyfile = '---------MEPhIST.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'MEPhIST'

#%% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()

levels = 20
# levels = list(numpy.linspace(-60e-5, 9e-5, 100))

#%% Domain and mesh definition
domain = geometry.rectangle_domain(area=[0.05, 0.55, -0.4, 0.4])
mephist_vessel = geometry.mephist_vessel()

domain.set_subdomain(1, mephist_vessel)

geometry.generate_mesh_in_domain(domain=domain, density=100)

# plot(geometry.mesh)
# fu.save_contour_plot(PATH, '')

markers = MeshFunction("size_t", geometry.mesh, geometry.mesh.topology().dim(), geometry.mesh.domains())

# plot(markers)
# fu.save_contour_plot(PATH, 'Markers')

#%% Define function space and step coefficients
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh

u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(1)))

[p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_pow_2(p_correction=1e2, F_correction=1)

u = TrialFunction(V) # u must be defined as function before expression def

a = dot(grad(u)/r, grad(r_2*v))*dx - (p_coeff*r*r + F_2_coeff)*u*r*v*dx(1)
L = sum(point_sources)*r*v*dx(0)

# a = dot(grad(u)/r, grad(r_2*v))*dx - (p_coeff*r*r + F_2_coeff)*u*r*v*dx(0)
# L = Constant(0)*r*v*dx 

u = Function(V)
solve(a == L, u, bc)

#%% Post solve
fu.What_time_is_it(t0, 'Variational problem solved')
fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')

fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
logger.log_n_output_colored_message(colored_message="'Done'\n", color='red', white_message='')
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
import numpy
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------Spheromak.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'Spheromak'

#%% Needed objects and contour levels
boundary = BoundaryConditions()
geometry = Geometry()
mesh_density = 50

levels = 20
# levels = numpy.linspace(-0.14, 0.15, 25)

#%% Domain and mesh definition
# geometry.rectangle_mesh_init(r1 = 0.05, r2 = 0.7, z1 = -0.6, z2 = 0.6, default_mesh = 80)
geometry.arbitrary_mesh_init()
fu.plot_mesh(geometry.mesh, PATH)

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
boundary.spheromak_boundary_condition(psi_0 = 0.1, R = 0.3, alpha = 1.1)
u_D = boundary.psi_sol_expr
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

a = dot(grad(u)/r, grad(r_2*v))*dx
L = boundary.spheromak_right_hand_expr * r * v * dx

u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

# fu.countour_plot_via_mesh(geometry, interpolate(boundary.psi_sol_expr, V), levels = levels, PATH = PATH, plot_title = '')
fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')
logger.info("'Done'"+"\n")
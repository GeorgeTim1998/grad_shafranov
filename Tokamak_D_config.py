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
current_pyfile = '---------Tokamak_D_config.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'Tokamak_D_config'

#%% Needed objects and contour levels
fu.What_time_is_it(t0, "Start problem")
boundary = BoundaryConditions()
geometry = Geometry()
mesh_density = 50

levels = 50
levels = numpy.linspace(-0.10, 10, 25)

#%% Domain and mesh definition
geometry.rectangle_mesh_init(r1 = 0.02, r2 = 1.5, z1 = -0.7, z2 = 0.7, default_mesh = 200)
fu.plot_mesh(geometry.mesh, PATH)
# geometry.arbitrary_mesh_init()

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
boundary.tokamak_D_config_boundary_condition(psi_0 = 0.1, eps = 0.05, a = 1, b = 0.5) #!!!
u_D = boundary.psi_sol_expr
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

a = dot(grad(u)/r, grad(r_2*v))*dx
L = boundary.tokamak_D_config_right_hand_expr * r * v * dx

u = Function(V)
solve(a == L, u, bc)

#%% Post solve
fu.What_time_is_it(t0, "Problem is solved")

# fu.countour_plot_via_mesh(geometry, interpolate(boundary.psi_sol_expr, V), levels = levels, PATH = PATH, plot_title = '')
fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

fu.ErrorEstimate(u = u, u_D = u_D, mesh = geometry.mesh)

logger.info("'Done'"+"\n")
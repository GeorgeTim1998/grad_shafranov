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
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------Spheromak.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'Spheromak'
#%% Needed objects
boundary = BoundaryConditions()
geometry = Geometry()
mesh_density = 50

#%% Domain and mesh definition
geometry.rectangle_mesh_init(r1 = 0.1, r2 = 0.5, z1 = -0.3, z2 = 0.3, default_mesh = 50)
fu.plot_mesh(geometry.mesh, PATH)
# domain = fu.spheromak_bounfdary(R, delta, 500)
# mesh = mshr.generate_mesh(domain, 50)

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
boundary.spheromak_boundary_condition(psi_0 = 1, R = 0.4, alpha = M.MEPhIST().delta)
u_D = boundary.psi_sol_expr
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

a = dot(grad(u)/r, grad(r_2*v))*dx
L = boundary.spheromak_right_hand_expr * r * v * dx

u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

fu.countour_plot_via_mesh(geometry, interpolate(boundary.psi_sol_expr, V), levels = 20, PATH = PATH, plot_title = '')
logger.info("'Done'"+"\n")
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
import numpy

#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------FRC.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'FRC'

#%% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()
levels = 20

psi_0 = 0.1
m0_p2 = 0.277
logger.log_n_output_colored_message(colored_message="m0 * p2 = ", color='green', white_message=str(m0_p2))
logger.log_n_output_colored_message(colored_message="psi_0 = ", color='green', white_message=str(psi_0))

#%% Domain and mesh definition
geometry.rectangle_mesh_init(r1 = 0.01, r2 = 1.04, z1 = -0.5, z2 = 0.5, default_mesh = 500) # Howell2014

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0.1")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

u = Function(V) 
a = dot(grad(u)/r, grad(r_2*v))*dx + 2*m0_p2 * r * r * u/psi_0 * r * v * dx
solve(a == 0, u, bc)

#%% Post solve
fu.What_time_is_it(t0, "Problem is solved")

fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

logger.info("'Done'"+"\n")

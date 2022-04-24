from fenics import *
import matplotlib.pyplot as matplt
from sympy import degree
import logger
import mshr
import time
import funcs as fu
import MEPHIST_data as M
import logger
from geometry import Geometry
from boundary_conditions import BoundaryConditions
import numpy
from expressions import Expressions

#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------1D_problem.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = '1D_problem'

#%% Needed objects and contour levels
fu.What_time_is_it(t0, "Start problem")
geometry = Geometry()
my_expressions = Expressions()

#%% All definition
geometry.interval_mesh_init(a=0.1, b=2, default_mesh=150)
my_expressions.axissymm_config(j0=1, r0=geometry.b)

#%% Function space and oundary conditione
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

boundary_value = Expression('0', degree=1)
logger.info("boundary value = %s" % boundary_value._cppcode)

bc = DirichletBC(V, boundary_value, "near(x[0], %f)" % geometry.b)

def boundary(x, on_boundary):
    return near(x[0], geometry.b)

#%% Problem definition
a = grad(u)[0] * v * dx
L = my_expressions.axissymm_config_right_hand_expr * v * dx

#%% Solvelog_n_output_colorlog_n_output_colored_messageed_message
u = Function(V)
solve(a == L, u, bc)

#%% Plot
fu.plot_1D(PATH, u, geometry)
fu.ErrorEstimate(u=u, u_D=interpolate(my_expressions.axissymm_config_solution_expr, V=V), mesh=geometry.mesh)

logger.info("'Done'" + "\n")
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
from OneD_problem_params import Problem

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
p = Problem()

#%% All definition
# for default_mesh in p.default_mesh_array:
geometry.interval_mesh_init(a=p.a, b=p.b, default_mesh=p.default_mesh)
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
fu.What_time_is_it(t0, "Problem solved")

#%% Plot
fu.plot_1D(PATH, u, geometry)
p.errors.append(fu.ErrorEstimate(u=u, u_D=interpolate(my_expressions.axissymm_config_solution_expr, V=V), mesh=geometry.mesh)[1])
# fu.plot_error_vs_mesh_density(p.default_mesh_array, p.errors, PATH)
# fu.save_errors_to_file(p.default_mesh_array, p.errors, PATH)

logger.info("'Done'" + "\n")
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
from FRC_problem_params import Problem

#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------FRC.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'FRC'

#%% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()
problem = Problem()

#%% Domain and mesh definition
geometry.rectangle_mesh_init(r1 = problem.domain_geometry[0], r2 = problem.domain_geometry[1], z1 = problem.domain_geometry[2], z2 = problem.domain_geometry[3], default_mesh = problem.default_mesh) # Howell2014

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition(problem.boundary_condition_str)
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

logger.log_n_output_colored_message(colored_message="m0 * p2 = ", color='green', white_message=str(problem.m0_p2))
logger.log_n_output_colored_message(colored_message="psi_0 = ", color='green', white_message=str(problem.psi_0))

u = Function(V) 
a = dot(grad(u)/r, grad(r_2*v))*dx - 2*problem.m0_p2 * r * r * u/problem.psi_0 * r * v * dx
solve(a == 0, u, bc)

#%% Post solve
fu.What_time_is_it(t0, "Problem is solved")

fu.countour_plot_via_mesh(geometry, u, levels = problem.contour_levels, PATH = PATH, plot_title = '')
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

logger.info("'Done'"+"\n")

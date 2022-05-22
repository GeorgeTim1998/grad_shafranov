### В этом файле решаю задачу о восстановлении равновевсия в FRC (Field Reverced Configuration)
#%% Imports
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
# geometry.rectangle_mesh_init(r1 = problem.domain_geometry[0], r2 = problem.domain_geometry[1], z1 = problem.domain_geometry[2], z2 = problem.domain_geometry[3], default_mesh = problem.default_mesh) # Howell2014

domain = geometry.circle_domain(centre_point=problem.centre_point, radius=problem.radius, segments=problem.segments)
geometry.generate_mesh_in_domain(domain, density=problem.mesh_density)

#%% Define function space and
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition(problem.boundary_condition_str)
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

u = Function(V) 
# u = interpolate(Expression("1e6", degree=1), V)
a = dot(grad(u)/r, grad(r_2*v))*dx - 2*problem.m0_p2 * r * r * u/problem.psi_0 * r*v*dx

solve(a == 0, u, bc)

#%% Post solve
fu.What_time_is_it(t0, "Problem is solved")

fu.countour_plot_via_mesh(geometry, u, levels = problem.contour_levels, PATH = PATH, plot_title = '')
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

logger.info("'Done'"+"\n")

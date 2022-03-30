#%% Imports
from matplotlib import interactive
from sympy import true
from funcs import print_colored
from imports import *
import time
import logger
import logging
import mshr
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------GS_MEPhIST_nomshr.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
logging.getLogger('FFC').setLevel(logging.WARNING)
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'MEPhIST'
#%% Geometry
default_mesh = fu.DEFAULT_MESH
abs_tol, rel_tol = 1e-10, 1e-9 # default = 1e-10, 1e-9
maximum_iterations = 100

logger.info('Newton Solver params: abs_tol = %e, rel_tol = %e, max_iter = %d' % (abs_tol, rel_tol, maximum_iterations))

eps = fu.EPS # when zero maybe inf (1/r)
r1, z1 = fu.R1 + eps, fu.Z1 # see Krat's unpublishet article
r2, z2 = fu.R2, fu.Z2
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]

logger.log_n_output('DEFAULT_MESH = %d' % default_mesh, 'green')
logger.info('(R1 +) EPS = %f' % eps)
logger.info('R1 = %f, Z1 = %f' % (r1, z1))
logger.info('R2 = %f, Z2 = %f' % (r2, z2))

mesh_r, mesh_z = default_mesh, abs(int(default_mesh * (z2-z1)/(r2-r1)))
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
logger.info("Number of cells: %d, Number of vertices: %d" % (mesh.num_cells(), mesh.num_vertices()))

V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
#%% Define funcs and weights
u = Function(V) # u must be defined as function before expression def
v = TestFunction(V)

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such
#%% Boundary conditions and function space V
u_D_str = '0'
u_D = Expression(u_D_str, degree = 1) # Define boundary condition
boundary = fu.DIRICHLET_BOUNDARY
logger.info('u_D = %s' % u_D_str)

if boundary == fu.DIRICHLET_BOUNDARY:
    bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле
    logger.info(fu.DIRICHLET_BOUNDARY)
else:
    bc = DirichletBC(V, u_D, fu.Neumann_boundary) #гран условие как в задаче дирихле
    logger.info(fu.NEUMANN_BOUNDARY)
#%% Problem pre-solve
todo = fu.SOLVE_PLASMA_POINT_SOURCES_EXPLICIT
logger.info("We are doing: %s" % str(todo))
fu.print_colored(todo, 'green')

p_pow = 2
F_pow = 2
logger.info("p_pow = %d, F_pow = %d" % (p_pow, F_pow))
f_text = fu.Hand_input(p_pow, F_pow)

f_expr = Expression(f_text, u = u, degree = 2)
u = fu.Initial_guess_for_u(u, 0)

A1 = 1 
A2 = 10
step = 1

alpha_array = numpy.linspace(A1, A2, 1+int((A2-A1)/step))
alpha_array = [A1]  

fu.What_time_is_it(t0, 'Problem posted')
#%% Problem solve
for alpha in alpha_array:
    point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(alpha)))
    L = sum(point_sources)*r*v*dx 

    if todo == fu.SOLVE_PLASMA_POINT_SOURCES:
        a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx - L 
        solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol, "maximum_iterations": maximum_iterations}})
    elif todo == fu.SOLVE_PLASMA:
        a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx
        solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol, "maximum_iterations": maximum_iterations}})
    else:
        if todo == fu.SOLVE_PLASMA_POINT_SOURCES:
            u = TrialFunction(V)
            a = dot(grad(u)/r, grad(r_2*v))*dx 
            u = Function(V)
            solve(a == L, u, bc)
        elif todo == fu.SOLVE_PLASMA_POINT_SOURCES_EXPLICIT:
            u = TrialFunction(V)
            a = dot(grad(u)/r, grad(r_2*v))*dx - (8.0*u*r_2 - 16.0*u)*r*v*dx
            u = Function(V)
            solve(a == L, u, bc)
            # a = dot(grad(u)/r, grad(r_2*v))*dx - (8.0*u*r_2 - 16.0*u)*r*v*dx - L 
            # solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol, "maximum_iterations": maximum_iterations}})
#%% Endproblem
    fu.What_time_is_it(t0, 'Variational problem solved')
    logger.log_n_output("(\u03C3*)\u03B1 = %e" % alpha, 'green')
    fu.print_colored("Solve for p_pow = %s, F_pow = %s" % (p_pow, F_pow), 'green')
    fu.print_colored("Calculations, however bad, finished", 'green')

    plot_title = "\u03C3*%.2e, #%d" % (alpha, int(default_mesh))
    fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, '', [mesh_r, mesh_z], plot_title, 20)
    fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
    logger.info("'Done'"+"\n")

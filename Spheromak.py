from fenics import *
import matplotlib.pyplot as matplt
import logger
import mshr
import time
import funcs as fu
import MEPHIST_data as M
import logger
from geometry import Geometry
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------Spheromak.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'Spheromak'
#%% Spheromak geometry and initial values
R = 0.3
logger.log_n_output("R = %f" % R, "green")

delta = M.MEPhIST().delta
logger.log_n_output("Вытянутость = %f" % delta, "green")

psi_0 = 1
logger.log_n_output("psi_0 = %f" % psi_0, "green")

mesh_density = 50
#%% Boundary conditions
u_D_str = '0'
u_D = Expression(u_D_str, degree = 1) # Define boundary condition
logger.info('u_D = %s' % u_D_str)
#%% Domain and mesh definition
geometry = Geometry()
geometry.rectangle_mesh_init(r1 = 0.1, r2 = 0.5, z1 = -0.3, z2 = 0.3, default_mesh = 100)
# domain = fu.spheromak_boundary(R, delta, 500)
# mesh = mshr.generate_mesh(domain, 50)
fu.plot_mesh(geometry.mesh, PATH)
#%% Formulate problem
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh
u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

f_expr = Expression("pow(x[0], 2) / pow(%s, 4) * (1 + pow(%s, 2)) * %s" % (R, delta, psi_0), degree = 2)
logger.log_n_output("Rigth hand part:", 'red')
logger.log_n_output(f_expr._cppcode, 'white')

[r_2, r] = geometry.operator_weights(V)

a = dot(grad(u)/r, grad(r_2*v))*dx
L = f_expr * r * v * dx
#%% Boundary conditions and solve
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")

fu.countour_plot_via_mesh(geometry, u, levels = 20, PATH = PATH, plot_title = '')
logger.info("'Done'"+"\n")
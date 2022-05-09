#%% imports
from fenics import *
import matplotlib.pyplot as plt
import datetime
import numpy 
import sympy
from termcolor import colored
import pylab as plt
import funcs as fu
from geometry import Geometry

PATH = 'Analytical'
CONTOUR_AMOUNT = 25
#%% Functions
def Form_f_text(A1, A2):
    #A1 = 4*pi*p', A2 = FF'
    #deriviation are calculated using sympy library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("INVERCED right-hand equation side: \n", 'magenta') + f_text)

    return f_text

def Analyt_sol(c, A1, A2):
    x = sympy.symbols('x[0]') # r coordinate
    z = sympy.symbols('x[1]') # r coordinate
    #sympy.log
    psi_p = A1 * pow(x, 4) + A2 * pow(z, 2) #private solution
    psi_gen = \
        c[0] + \
        c[1] * pow(x, 2) + \
        c[2] * (pow(x, 4) - 4*pow(x, 2)*pow(z, 2)) + \
        c[3] * (-pow(z, 2)) # general solution the rest of the 4th term is defined in MyLog(c) func
        #c[3] * (pow(x, 2)*sympy.log(x)- pow(z, 2)) # general solution
        #pow(x, 2)*sympy.log(x) 
    
    my_log = MyLog(c)
    
    psi_text = sympy.printing.ccode(psi_p + psi_gen)
    psi_p_text = sympy.printing.ccode(psi_p)
    
    # final_sol = psi_text 
    final_sol = psi_text + ' + ' + my_log
    print(colored("Private solution: \n", 'magenta') + psi_p_text)
    print(colored("Analytical solution: \n", 'magenta') + final_sol)
    # print(colored("Analytical solution: \n", 'magenta') + psi_text)

    return final_sol

def MyLog(c):
    x = sympy.symbols('x[0]') # r coordinate
    pre_log = c[3] * pow(x, 2)
    
    pre_log_text = sympy.printing.ccode(pre_log)
    log_text = "%s*std::log(%s)" % (pre_log_text, 'x[0]') # assemble function of the point source
    print(colored("Problem term in analyt solution: \n", 'magenta') + log_text)
    #c[3] * (pow(x, 2)*log(x)) # general solution
    
    return log_text

def ErrorEstimate(u, u_D, mesh):
    # Compute error in L2 norm
    error_L2 = errornorm(u_D, u, 'L2')

    # Compute maximum error at vertices
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = numpy.max(numpy.abs(vertex_values_u_D - vertex_values_u))

    # Print errors
    print(colored('error_L2  = ', 'red'), error_L2)
    print(colored('error_max = ', 'red'), error_max)
    
    return error_L2, error_max

print(colored("GS_Soloviev_analyt.py", 'green'))

#%% Needed objects
geom = Geometry()

#%% paremeters definition
area = [0.2, 2.2, -1, 1] # format is: [r1, r2, z1, z2]
geom.rectangle_mesh_init(area[0], area[1], area[2], area[3], default_mesh=150)

A1, A2 = 0.14, -0.01 # values from Ilgisonis2016, 244
c = [1, -0.22, -0.01, -0.08] # values from Ilgisonis2016, 244

f_text = Form_f_text(-8 * A1, -2 * A2) # form right hand side that corresponds to analytical solution
psi_text = Analyt_sol(c, A1, A2) # см. научка.txt. Там есть вывод, как и куда надо подставлять


#%% Create mesh and define function space
V = FunctionSpace(geom.mesh, 'P', 1) # standard triangular mesh
u_D = Expression(psi_text, degree = 4) # Define boundary condition

psi = interpolate(u_D, V) #plot exact solution
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле
#%% Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f_expr = Expression(f_text, degree = 2)

[r_2, r] = geom.operator_weights(V)

a = dot(grad(u)/r, grad(r_2*v))*dx
L = f_expr*r*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)
#%% Compute errors
[err_L2, err_max] = ErrorEstimate(u, psi, geom.mesh)
# fu.Write2file_errors(mesh_r, mesh_z, err_L2, err_max)
#%% Save output
levels = 40
fu.countour_plot_via_mesh(geom, u, levels, PATH=PATH, plot_title='')
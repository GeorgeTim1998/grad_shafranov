#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, isinteractive, title
import datetime
import numpy 
import sympy
from termcolor import colored
import pylab
import funcs as fu
#%% Functions
def Save_figure(addition):
    # Plot solution and mesh. Save plot
    #nothing passed to function, because variables are global
    if plot_mesh == 1:
        plot(mesh)

    mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
    
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  #get current time to make figure name unique

    #create a path to save my figure to. For some reason now I cant save using relative path
    path_my_file = '/home/george/Projects2/Projects/Figures/Bourder/' + time_title

    if addition == '_notitle':
        plt.savefig(path_my_file + addition + '.png', dpi = dpi) #no title figure for reports
    elif addition == '_title':
        plt.title('Analyt Soloviev: ' + mesh_title + "\n" + f_expr._cppcode) # titled figure for my self
        plt.savefig(path_my_file + addition + '.png', dpi = dpi)
    else:
        plt.title(addition) # titled figure for my self
        plt.savefig(path_my_file + addition + '.png', dpi = dpi)

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

print(colored("GS_Soloviev_new.py", 'green'))
#%% paremeters definition
r0, z0 = 100, 0 # starting point for calculations
square = 1 # square size
mesh_r, mesh_z = 200*square, 200*square # mesh for r-z space
r1, z1 = r0 - 0.5*square, z0 - 0.5*square
r2, z2 = r0 + 0.5*square, z0 + 0.5*square
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]

rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

del r1, r2, z1, z2

plot_mesh = 0 #choose whether to plot mesh or not
save_NoTitle = 0 #save figure that doesnt have title in it
show_plot = 0 # show plot by the end of the program or not
dpi = 200 # quality of a figure 

A1, A2 = 0.14, 0.01 # values from Ilgisonis2016, 244
f_text = fu.Form_f_text(A1, A2) # form right hand side that corresponds to analytical solution
#%% Create mesh and define function space
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 1) # Define boundary condition

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле
#%% Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f_expr = Expression(f_text, degree = 2)
r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot(grad(u)/r, grad(r_2*v))*dx
L = f_expr*r*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)
fig = plot(u) # its fenics' plot not python's
pylab.colorbar(fig)
#%% Save output
Save_figure('_title')
vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
vtkfile << u
#%% 'plt.show()' holds plot while the programm is still running
if show_plot == 1:
    plt.show()
#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
from numpy import add
import sympy
from termcolor import colored
#%% Functions
def Form_f_text(A1, A2):
    #A1 = 4*pi*p', A2 = FF'
    #deriviation are calculated using sympy library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("\nINVERCED right-hand equation side: ", 'magenta') + f_text)

    return f_text

def Save_figure(addition):
    # Plot solution and mesh. Save plot
    #nothing passed to function, because variables are global
    if plot_mesh == 1:
        plot(mesh)

    mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
    
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  #get current time to make figure name unique

    #create a path to save my figure to. For some reason now I cant save using relative path
    path_my_file = '/home/george/Projects2/Projects/Figures/Analytical/' + time_title

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
        c[3] * (pow(x, 2)*sympy.log(x)- pow(z, 2)) # general solution
        #pow(x, 2)*sympy.log(x) 
    psi_text = sympy.printing.ccode(psi_p + psi_gen)
    psi_p_text = sympy.printing.ccode(psi_p)
    print(colored("\nAnalytical solution: ", 'magenta') + psi_text + "\n")
    print(colored("Private solution: ", 'magenta') + psi_p_text + "\n")

    return psi_text
#%% paremeters definition
mesh_r, mesh_z = 100, 100 # mesh for r-z space
area = [0.2, 2.2, -1, 1] # format is: [r1, r2, z1, z2]
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

plot_mesh = 0 #choose whether to plot mesh or not
save_NoTitle = 0 #save figure that doesnt have title in it
show_plot = 0 # show plot by the end of the program or not
dpi = 200 # quality of a figure 


A1, A2 = 0.14, -0.01
c = [1, -0.22, -0.01, -0.08] #coefficients used for analytical solution

f_text = Form_f_text(A1 * 8, A2 * 2) # form right hand side that corresponds to analytical solution
psi_text = Analyt_sol(c, -A1, -A2) #A1&A2 defined as: A1 = 4*pi*p', A2 = FF'!!!
#%% Create mesh and define function space
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression(psi_text, degree = 3) # Define boundary condition

psi = interpolate(u_D, V) #plot exact solution
plot(psi)
Save_figure('_Analit')

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле
#%% Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f_expr = Expression(f_text, degree = 2)
w = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot(grad(u)/w, grad(w*v))*dx
L = f_expr*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)
plot(u) # its fenics' plot not python's
#%% Save output
Save_figure('_title')
vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
vtkfile << u
#%% 'plt.show()' holds plot while the programm is still running
if show_plot == 1:
    plt.show()
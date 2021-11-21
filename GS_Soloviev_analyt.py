#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
#%% Functions
def Form_f_text(A1, A2):
    #A1 = -4*pi*p', A2 = -FF'
    #expression is inverced because f deined as -f0 (what you see in GS equation)
    #deriviation are calculated using sympy library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print("\n" + f_text + "\n")

    return f_text

def Save_figure(mesh, mesh_r, mesh_z, plot_mesh, save_NoTitle, dpi):
    # Plot solution and mesh. Save plot

    if plot_mesh == 1:
        plot(mesh)

    mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
    
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  #get current time to make figure name unique

    #create a path to save my figure to. For some reason now I cant save using relative path
    path_my_file = '/home/george/Projects2/Projects/Figures/' + time_title

    if save_NoTitle != 0:
        plt.savefig(path_my_file + '_notitle.png', dpi = dpi) #no title figure for reports

    plt.title('Soloviev: ' + mesh_title + "\n" + f_expr._cppcode) # titled figure for my self
    plt.savefig(path_my_file + '_title.png', dpi = dpi)
#%% paremeters definition
mesh_r, mesh_z = 100, 100 # mesh for r-z space
rect_low = Point(0, -1) #define rectangle size
rect_high = Point(1, 1)

plot_mesh = 0 #choose whether to plot mesh or not
save_NoTitle = 0 #save figure that doesnt have title in it
dpi = 200 # saved figures quality
show_plot = 0

A1, A2 = 1, 2
f_text = Form_f_text(A1, A2)
#%% Create mesh and define function space
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 0) # Define boundary condition #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле
#%% Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f_expr = Expression(f_text, degree = 2)
f1 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # comment in {}

a = dot(grad(u)/f1, grad(f1*v))*dx
L = f_expr*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)
plot(u) # its fenics' plot not python's
#%% Save output
Save_figure(mesh, mesh_r, mesh_z, plot_mesh, save_NoTitle, dpi)
vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
vtkfile << u
#%% 'plt.show()' holds plot while the programm is still running
if show_plot == 1:
    plt.show()
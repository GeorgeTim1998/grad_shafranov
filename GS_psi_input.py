#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
from termcolor import colored
#%% Functions
def Save_figure(addition, f_expr):
    # Plot solution and mesh. Save plot
    #nothing passed to function, because variables are global
    if plot_mesh == 1:
        plot(mesh)

    mesh_title = "%sx%s mesh" % (str(mesh_r), str(mesh_z))
    curr_time = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(curr_time)  #get current time to make figure name unique
    path_my_file = '/home/george/Projects2/Projects/Figures/Psi_input/%s' % time_title


    if addition == '_notitle':
        plt.savefig("%s%s.png" %(path_my_file, addition), dpi = dpi) #no title figure for reports
    elif addition == '_title':
        plt.title('Point source: %s\n%s' % (mesh_title, f_expr._cppcode)) # titled figure for my self
        plt.savefig("%s%s.png" %(path_my_file, addition), dpi = dpi) #no title figure for reports
    else:
        plt.title(addition) # titled figure for my self
        plt.savefig("%s%s.png" %(path_my_file, addition), dpi = dpi) #no title figure for reports
#%% paremeters, mesh, V, u definition
mesh_r, mesh_z = 100, 100 # mesh for r-z space
area = [0, 1, -1, 1] # format is: [r1, r2, z1, z2]
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

plot_mesh = 0 #choose whether to plot mesh or not
save_NoTitle = 0 #save figure that doesnt have title in it
dpi = 200 # saved figures quality

mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u = Function(V) # for non linear equations 'u' must be defined via Function() 
#%% Equation's right hand definition 
# expression is inverced because f deined as -f0 (what you see in GS equation)
#!!!calc deriviation using sympy
psi = sympy.symbols('u') # flux function #think tomorrow how to define argument psi!
x = sympy.symbols('x[0]') # r coordinate. used for easy writing of expressions

p_psi = pow(psi, 2) #pressure function
F_psi = pow(psi, 1) # poloidal current function

dp_psi = sympy.diff(p_psi, psi) #pressure and F deriviation
dF_psi = sympy.diff(F_psi, psi) #compiler breaks when 

f_text = 4 * pi * pow(x, 2) * dp_psi + F_psi*dF_psi #right hand expression
#f_text = 4 * pi * pow(x, 2) * psi + 1 #right hand expression
f_text = sympy.printing.ccode(f_text)
print(colored("Right-hand side: ", 'magenta') + f_text)
f_expr = Expression(f_text, u = u, degree = 2)
#%% Create problem
#  mesh, u and V defined above

# Define boundary condition
u_D = Expression('0', degree=0) #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

# Define variational problem
v = TestFunction(V)
w = interpolate(Expression('x[0]*x[0]', degree = 2), V) # comment in {}

a = dot( grad(u)/w, grad(w*v) )*dx
L = f_expr*v*dx

# Compute solution
solve(a - L == 0, u, bc)
#%% create a path to save my figure to. For some reason now I cant save using relative path
mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
time_title = str(ttime)  #get current time to make figure name unique

plot(u, tol = 1e-6 ) # its fenics' plot not python's
Save_figure('_title', f_expr)

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u
#%% hold plot to show. Programm is still running
#plt.show()
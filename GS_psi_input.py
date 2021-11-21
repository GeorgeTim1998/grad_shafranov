#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
#%% paremeters, mesh, V, u definition
mesh_r = 10 #mesh
mesh_z = 10
plot_mesh = 0 #choose whether to plot mesh or not
rectangle_low = Point(0, -1) #define rectangle size
rectangle_high = Point(1, 1)

mesh = RectangleMesh(rectangle_low, rectangle_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u = Function(V) # for non linear equations 'u' must be defined via Function() 

save_NoTitle = 0 #save figure that doesnt have title in it
dpi = 200 # saved figures quality
#%% Equation's right hand definition 
# expression is inverced because f deined as -f0 (what you see in GS equation)
#!!!calc deriviation using sympy
psi = sympy.symbols('u') # flux function #think tomorrow how to define argument psi!
x = sympy.symbols('x[0]') # r coordinate. used for easy writing of expressions

p_psi = pow(psi, 2) #pressure function
F_psi = pow(psi, 2) # poloidal current function

dp_psi = sympy.diff(p_psi, psi) #pressure and F deriviation
dF_psi = sympy.diff(F_psi, psi) #compiler breaks when 

f_text = 4 * pi * pow(x, 2) * dp_psi + F_psi*dF_psi #right hand expression
#f_text = 4 * pi * pow(x, 2) * psi + 1 #right hand expression
f_text = sympy.printing.ccode(f_text)
print("\n" + f_text + "\n")
f = Expression(f_text, u = u, degree = 2)
#%% Create problem
#  mesh, u and V defined above

# Define boundary condition
u_D = Expression('0', degree=0) #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

# Define variational problem
v = TestFunction(V)
r2weigth = interpolate(Expression('x[0]*x[0]', degree = 2), V) # comment in {}

a = dot( grad(u)/r2weigth, grad(r2weigth*v) )*dx
L = f*v*dx

# Compute solution
solve(a - L == 0, u, bc)
#%% create a path to save my figure to. For some reason now I cant save using relative path
mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
time_title = str(ttime)  #get current time to make figure name unique

plot(u) # its fenics' plot not python's
if plot_mesh == 1:
    plot(mesh)

path_my_file = '/home/george/Projects2/Projects/Figures/' + time_title

if save_NoTitle != 0:
    plt.savefig(path_my_file + '_notitle.png', dpi = dpi) #no title figure for reports
plt.title('Psi_input: ' + mesh_title + "\n" + f._cppcode) # titled figure for my self
plt.savefig(path_my_file + '_title.png', dpi = dpi)

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u
#%% hold plot to show. Programm is still running
#plt.show()
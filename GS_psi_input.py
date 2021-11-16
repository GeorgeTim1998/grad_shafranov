#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
#from ufl.constantvalue import ConstantValue
#%% paremeters definition
mesh_r = 100 #mesh
mesh_z = 100
plot_mesh = 0 #choose whether to plot mesh or not
rectangle_low = Point(0, -1) #define rectangle size
rectangle_high = Point(1, 1)

save_NoTitle = 0 #save figure that doesnt have title in it
dpi = 200 # saved figures quality
#%% Equation's right hand definition 
# expression is inverced because f deined as -f0 (what you see in GS equation)
#!!!calc deriviation using sympy
psi = sympy.symbols('u') # flux function #think tomorrow how to define argument psi!
x = sympy.symbols('x[0]') # r coordinate. used for easy writing of expressions

p_psi = pow(psi, 0) #pressure function
F_psi = pow(psi, 1) # poloidal current function

dp_psi = sympy.diff(p_psi, psi) #pressure and F deriviation
dF_psi = sympy.diff(F_psi, psi) #compiler breaks when 
#u is present on the right side of equation.
# посмотри про переменную degree!!!

f_text = 4 * pi * pow(x, 2) * dp_psi + F_psi * dF_psi #right hand expression
f_text = sympy.printing.ccode(f_text)
print("\n" + f_text + "\n")

mesh = RectangleMesh(rectangle_low, rectangle_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
#u = TrialFunction(V) 
u = Function(V) # for non linear equations 'u' must be defined via Function() 
f = Expression(f_text, u = u, degree = 2)

#f_text = sympy.simplify(f_text) #make expression simpler then it is
#x2 = 2*pi*pow(x, 3)+x*10+1/(x+1) #examples!
#!!!calculate deriviation using sympy
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
{
    #r2weigth is r^2 that appears during
    #this function is used to define operator just like in G-Sh equation
    #make r2weigth available for spacial deriviations
    #if I want to state the problem as in
    # my "Вывод ур-ия Г-Ш.pdf" file then I probably need 'Expression' func
    # just like Ineed it to define right hand side of the equation
}

a = dot( grad(u)/r2weigth, grad(r2weigth*v) )*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a - L == 0, u, bc)
#%% create a path to save my figure to. For some reason now I cant save using relative path
# Plot solution and mesh. Save plot
mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
time_title = str(ttime)  #get current time to make figure name unique

plot(u) # its fenics' plot not python's
if plot_mesh == 1:
    plot(mesh)

path_my_file = '/home/george/Projects2/Projects/Figures/' + time_title

if save_NoTitle != 0:
    plt.savefig(path_my_file + '_notitle.png', dpi = dpi) #no title figure for reports
plt.title('Soloviev: ' + mesh_title + "\n" + f._cppcode) # titled figure for my self
plt.savefig(path_my_file + '_title.png', dpi = dpi)

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u
#%% hold plot to show. Programm is still running
#plt.show()
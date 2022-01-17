#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
#%% paremeters definition
mesh_r = 100 #mesh
mesh_z = 100
plot_mesh = 0 #choose whether to plot mesh or not
rectangle_low = Point(0, -1) #define rectangle size
rectangle_high = Point(1, 1)

save_NoTitle = 1 #save figure that doesnt have title in it
dpi = 200 # saved figures quality

# expression is inverced because f deined as -f0 (what you see in GS equation)
#!!!calc deriviation using sympy
A1 = 1#-4*pi*p'
A2 = 2#-FF'
x = sympy.symbols('x[0]') # r coordinate
f_1 = A1 * pow(x, 2) + A2
f_1 = sympy.printing.ccode(f_1)

print("\n" + f_1 + "\n")
#%% Create mesh and define function space
mesh = RectangleMesh(rectangle_low, rectangle_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh

# Define boundary condition
u_D = Expression('0', degree=0) #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Expression(str(A1) + '*x[0]*x[0]+' + str(A2), degree = 2)
f = Expression(f_1, degree = 2)
f1 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # comment in {}

a = dot(grad(u)/f1, grad(f1*v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh. Save plot
mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
A1_title = 'A1 = ' + str(A1) 
A2_title = 'A2 = ' + str(A2)
ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
time_title = str(ttime)  #get current time to make figure name unique


plot(u) # its fenics' plot not python's
if plot_mesh == 1:
    plot(mesh)


#create a path to save my figure to. For some reason now I cant save using relative path
path_my_file = '/home/george/Projects2/Projects/Figures/' + time_title

# if save_NoTitle != 0:
#     plt.savefig(path_my_file + '_notitle.png', dpi = dpi) #no title figure for reports
# plt.title('Soloviev: ' + mesh_title + "\n" + f._cppcode) # titled figure for my self
# plt.savefig(path_my_file + '_title.png', dpi = dpi)

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u
#%% hold plot to show. Programm is still running
#plt.show()
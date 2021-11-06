#%% main code
## in terminal do 'conda start phenics project'
# comment lines: cntr+K+C
# uncomment lines: cntr+K+U

from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
import datetime

#paremeters definition
mesh_r = 10 #mesh
mesh_z = 10
A1 = 1 #4*pi*p'
A2 = 1 #FF'

# Create mesh and define function space
mesh = UnitSquareMesh(mesh_r, mesh_z)
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh

# Define boundary condition
u_D = Expression('0', degree=0) #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #задача дирихле

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression(str(A1) + '*x[0]*x[0]+' + str(A2), degree = 2) 
#f = Expression('x[0]*x[0]+1', degree = 2) 
f1 = interpolate(Expression('x[0]*x[0]', degree = 2), V) #this function is used to define operator just like in G-Sh equation
#make f1 available for spacial deriviations
#if I want to state the problem as in
# my "Вывод ур-ия Г-Ш.pdf" file then I probably need 'Expression' func
# just like Ineed it to define right hand side of the equation

a = dot(grad(u)/f1, grad(f1*v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh. Save plot
mesh_title = str(mesh_r) + 'x' + str(mesh_z) + ' mesh'
A1_title = 'A1 = ' + str(A1) 
A2_title = 'A2 = ' + str(A2)
ttime = datetime.datetime.now().strftime("%d%m%Y:%H%M%S")
time_title = str(ttime)  #get current time to make figure name unique

plot(u)
plot(mesh)

#plt.savefig('Figures/' + time_title + '_no_title') #no title figure for reports
plt.title(mesh_title + "\n" + A1_title + ', ' + A2_title) # titled figure for my self
plt.savefig('Figures/' + time_title + '_title')

# Save solution to file in VTK format
vtkfile = File('poisson/solution.pvd')
vtkfile << u
#%% compare with exact solution
  
# # Compute error in L2 norm
# error_L2 = errornorm(u_D, u, 'L2')

# # Compute maximum error at vertices
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)
# import numpy as np
# error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# # Print errors
# print('error_L2  =', error_L2)
# print('error_max =', error_max)

#%% hold plot to show. Programm is still running
#plt.show()
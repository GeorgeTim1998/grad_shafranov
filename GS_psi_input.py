#%% imports
from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, title
import datetime
import sympy
import funcs as fu
from termcolor import colored
import time

PATH = 'Psi_input'
t0 = time.time()
fu.What_time_is_it(t0, 'Start')
#%% paremeters, mesh, V, u definition
area = [0.4, 1, -1, 1] # format is: [r1, r2, z1, z2]
mesh_r, mesh_z = fu.DEFAULT_MESH, fu.DEFAULT_MESH # mesh for r-z space
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size [0, -1]x[1, 1]

V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u = Function(V) # for non linear equations 'u' must be defined via Function() 

psi = sympy.symbols('u') # flux function #think tomorrow how to define argument psi!
x = sympy.symbols('x[0]') # r coordinate. used for easy writing of expressions

p_psi = pow(psi, 2) #pressure function
F_psi = 1 - pow(psi, 2) # poloidal current function

dp_psi = sympy.diff(p_psi, psi) #pressure and F deriviation
dF_psi = sympy.diff(F_psi, psi) #compiler breaks when 

f_text = 4 * pi * pow(x, 2) * dp_psi + F_psi*dF_psi #right hand expression

p_equat_text = 4 * pi * pow(x, 2) * dp_psi
F_equat_text = F_psi * dF_psi

f_text = sympy.printing.ccode(f_text)
p_equat_text = sympy.printing.ccode(p_equat_text)
F_equat_text = sympy.printing.ccode(F_equat_text)


print(colored("Right-hand side: ", 'magenta') + f_text)
f_expr = Expression(f_text, u = u, degree = 2)
p_equat = Expression(p_equat_text, u = u, degree = 2)
F_equat = Expression(F_equat_text, u = u, degree = 2)

# Define boundary condition
u_D = Expression('0', degree=0) #psi flux is zero if you go far away

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

# Define variational problem
v = TestFunction(V)
r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot( grad(u)/r, grad(r_2*v) )*dx 
L = (f_expr)*v*dx
solve(a - L == 0, u, bc)

fu.What_time_is_it(t0, 'Plotting...')
fu.Contour_plot([area[0], area[1]], [area[2], area[3]], u, PATH, '', [mesh_r, mesh_z], '', 20)
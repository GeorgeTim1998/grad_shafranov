from imports import *
import time
#%% Pre-programm stuff
t0 = time.time()
print(colored("\n---------GS_Soloviev_new.py---------", 'green'))
print(colored("Date_Time is: %s" % fu.Time_name(), 'cyan'))
PATH = 'MEPhIST'
#%% Programm body

A1, A2 = 0.14, 0.01 # values from Ilgisonis2016, 244
f_text = fu.Form_f_text(A1, A2) # form right hand side that corresponds to analytical solution

mesh = 0
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 1) # Define boundary condition

def boundary(x, on_boundary):
    return on_boundary

fu.What_time_is_it(t0, 'Variational problem solved')

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

u = TrialFunction(V)
v = TestFunction(V)

f_expr = Expression(f_text, degree = 2)
point_sources = Expression(fu.ArrayOfPointSources(psd.PointSource()), degree = 2)

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot(grad(u)/r, grad(r_2*v))*dx
L = (sum(point_sources))*r*v*dx

print(colored("Default mesh = %d\n" % (mesh_size), 'green'))

u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, 'Variational problem solved')
# %% 
# fu.Write2file_umax_vs_def_mesh(mesh_r, mesh_z, fu.Twod_plot(u, r0, z1, z2, PATH))

fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")
# vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
# vtkfile << u
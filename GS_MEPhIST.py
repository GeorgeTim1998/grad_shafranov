from imports import *
import time
#%% Pre-programm stuff
t0 = time.time()
print(colored("\n---------GS_Soloviev_new.py---------", 'green'))
print(colored("Date_Time is: %s" % fu.Time_name(), 'cyan'))
PATH = 'MEPhIST'
plot_title = 'MEPhIST'
#%% Programm body
default_mesh = fu.DEFAULT_MESH
r1, z1 = 0, -0.4 # see Krat's unpublishet article
r2, z2 = 0.6, 0.4
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]

mesh_r, mesh_z = default_mesh, abs(int(default_mesh * (z2-z1)/(r2-r1)))

rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

A1, A2 = 0.14, 0.01 # values from Ilgisonis2016, 244
f_text = fu.Form_f_text(A1, A2) # form right hand side that corresponds to analytical solution

mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 1) # Define boundary condition

def boundary(x, on_boundary):
    return on_boundary

fu.What_time_is_it(t0, 'Variational problem solved')

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле

u = TrialFunction(V)
v = TestFunction(V)

f_expr = Expression(f_text, degree = 2)
point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource()))

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot(grad(u)/r, grad(r_2*v))*dx
L = sum(point_sources)*r*v*dx

print(colored("Default mesh = %d\n" % (default_mesh), 'green'))

u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, 'Variational problem solved')
# %% 
# fu.Write2file_umax_vs_def_mesh(mesh_r, mesh_z, fu.Twod_plot(u, r0, z1, z2, PATH))
fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, f_expr, [mesh_r, mesh_z], plot_title)
fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")
# vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
# vtkfile << u
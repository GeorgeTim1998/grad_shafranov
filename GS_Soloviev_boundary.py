from imports import *
import time
#%% Pre-programm stuff
t0 = time.time()
print(colored("\n---------GS_Soloviev_new.py---------", 'green'))
PATH = 'Border'
#%% paremeters definition
r0, z0 = 100, 0 # starting point for calculations
square = 0.25 # square size
default_mesh = 200
mesh_r, mesh_z = int(default_mesh*square), int(default_mesh*square) # mesh for r-z space
r1, z1 = r0 - 0.5*square, z0 - 0.5*square
r2, z2 = r0 + 0.5*square, z0 + 0.5*square
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]

rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

plot_mesh = 0 #choose whether to plot mesh or not
save_NoTitle = 0 #save figure that doesnt have title in it
show_plot = 0 # show plot by the end of the program or not

A1, A2 = 0.14, 0.01 # values from Ilgisonis2016, 244
f_text = fu.Form_f_text(A1, A2) # form right hand side that corresponds to analytical solution
#%% Create mesh and define function space
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 1) # Define boundary condition

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary) #гран условие как в задаче дирихле
#%% Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f_expr = Expression(f_text, degree = 2)
point_sources = Expression(fu.ArrayOfPointSources(psd.PointSource()), degree = 2)

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

a = dot(grad(u)/r, grad(r_2*v))*dx
L = (sum(point_sources))*r*v*dx
#%% Compute solution
print(colored("Default mesh = %f\nSquare size = %f" % (default_mesh, square), 'blue'))
u = Function(V)
solve(a == L, u, bc)
fu.What_time_is_it(t0, 'Variational problem solved')

# fu.Write2file_umax_vs_def_mesh(mesh_r, mesh_z, fu.Twod_plot(u, r0, z1, z2, PATH))
fu.Write2file_umax_vs_square_size(mesh_r, mesh_z, fu.Twod_plot(u, r0, z1, z2, PATH))
fu.What_time_is_it(t0, "Cross section plotted through r0 = %s" % r0)

fig = plot(u) # its fenics' plot not python's
pylab.colorbar(fig)
#%% Save output
fu.Save_figure(f_expr, mesh_r, mesh_z, '_title', PATH)
# vtkfile = File('poisson/solution.pvd') # Save solution to file in VTK format
# vtkfile << u

fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")
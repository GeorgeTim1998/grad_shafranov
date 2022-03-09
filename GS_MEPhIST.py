from imports import *
import time
#%% Pre-programm stuff
t0 = time.time()
print(colored("\n---------GS_MEPhIST.py---------", 'green'))
print(colored("Date_Time is: %s" % fu.Time_name(), 'cyan'))
PATH = 'MEPhIST'
plot_title = 'MEPhIST'
#%% Geometry
default_mesh = fu.DEFAULT_MESH
abs_tol, rel_tol = 1e-10, 8e-6 # default = 1e-10, 1e-9
eps = fu.EPS # when zero maybe inf (1/r)
r1, z1 = fu.R1 + eps, fu.Z1 # see Krat's unpublishet article
r2, z2 = fu.R2, fu.Z2
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]
mesh_r, mesh_z = default_mesh, abs(int(default_mesh * (z2-z1)/(r2-r1)))
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point

mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high
V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
#%% Boundary conditions and function space V
u_D_str = '1e-1'
u_D = Expression(u_D_str, degree = 1) # Define boundary condition
bc = DirichletBC(V, u_D, fu.Neumann_boundary) #гран условие как в задаче дирихле

u = Function(V) # u must be defined as function before expression def
v = TestFunction(V)

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such
point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource()))

fu.What_time_is_it(t0, 'Variational problem solved')
print(colored("Default mesh = %d\n" % (default_mesh), 'green'))
#%% SOLVING PROBLEM#1
# p_pow = 1
# F_pow = 2
# f_text = fu.Hand_input(p_pow, F_pow)

# f_expr = Expression(f_text, u = u, degree = 2)
# a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx
# solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}})
# fu.What_time_is_it(t0, "Solve for p_pow = %s, F_pow = %s" % (p_pow, F_pow))
# fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, '', [mesh_r, mesh_z], '', 20)

# fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, '', [mesh_r, mesh_z], '', 20)
# fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")
#%% SOLVING PROBLEM#2
p_pow = 2
F_pow = 2
f_text = fu.Hand_input(p_pow, F_pow)
u = fu.Initial_guess_for_u(u, 0)

f_expr = Expression(f_text, u = u, degree = 2)

L = sum(point_sources)*r*v*dx 
a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx - L 
solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}})

# a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx
# solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}})

fu.What_time_is_it(t0, "Solve for p_pow = %s, F_pow = %s" % (p_pow, F_pow))
#%% Endproblem
print(colored("Calculations, however bad, finished", 'green'))

fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, '', [mesh_r, mesh_z], u_D_str, 20)
fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")
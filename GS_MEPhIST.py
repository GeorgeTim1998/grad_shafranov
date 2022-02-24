from imports import *
import time
#%% Pre-programm stuff
t0 = time.time()
print(colored("\n---------GS_MEPhIST.py---------", 'green'))
print(colored("Date_Time is: %s" % fu.Time_name(), 'cyan'))
PATH = 'MEPhIST'
plot_title = 'MEPhIST'
#%% Programm body
default_mesh = fu.DEFAULT_MESH
abs_tol = 1e-10 # default = 1e-10
rel_tol = 1e-9 # default = 1e-9
eps = 0.005 # when zero maybe inf (1/r)
r1, z1 = 0 + eps, -0.4 # see Krat's unpublishet article
r2, z2 = 0.6, 0.4
area = [r1, r2, z1, z2] # format is: [r1, r2, z1, z2]
mesh_r, mesh_z = default_mesh, abs(int(default_mesh * (z2-z1)/(r2-r1)))
rect_low = Point(area[0], area[2]) #define rectangle size: lower point
rect_high = Point(area[1], area[3]) #define rectangle size: upper point
mesh = RectangleMesh(rect_low, rect_high, mesh_r, mesh_z) # points define domain size rect_low x rect_high

V = FunctionSpace(mesh, 'P', 1) # standard triangular mesh
u_D = Expression('0', degree = 1) # Define boundary condition

# def Dirichlet_boundary(x, on_boundary):
#     tol = 1e-10
#     return on_boundary and (abs(x[1] - z1) < tol or abs(x[1] - z2) < tol or abs(x[0] - r2) < tol)

def Dirichlet_boundary(x, on_boundary):
    return on_boundary

# class DirichletBoundary(SubDomain):
#     def inside(self, x, on_boundary):
#         tol = 1E-10   # tolerance for coordinate comparisons
#         return on_boundary and \
#                (abs(x[1] - z1) < tol or abs(x[1] - z2) < tol or abs(x[0] - r2) < tol)

# u0_boundary = DirichletBoundary()

fu.What_time_is_it(t0, 'Variational problem solved')

bc = DirichletBC(V, u_D, Dirichlet_boundary) #гран условие как в задаче дирихле

u = Function(V) # u must be defined as function before expression def
v = TestFunction(V)

f_text = fu.Hand_input()
f_expr = Expression(f_text, u = u, degree = 2)
point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource()))

r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such

print(colored("Default mesh = %d\n" % (default_mesh), 'green'))
# u = TrialFunction(V)

a = dot(grad(u)/r, grad(r_2*v))*dx - f_expr*r*v*dx

# a = dot(grad(u)/r, grad(r_2*v))*dx
# L = f_expr*r*v*dx

# L = sum(point_sources)*r*v*dx 

# u = Function(V)
#%% solve
# solve(a == 0, u, bc, "lu")
solve(a == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}})
# solve(a == 0, u, bc)
# solve(a - L == 0, u, bc)

print(colored("Calculations, however bad, finished", 'green'))

fu.Contour_plot([r1, r2], [z1,  z2], u, PATH, '', [mesh_r, mesh_z], '', 20)
fu.What_time_is_it(t0, "3D plot of \u03C8(r, z) is plotted")

# problem = NonlinearVariationalProblem(a == 0, u, bc)
# solver = NonlinearVariationalSolver(problem)
# solver_parameters={"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}
# solver.parameters.update(solver_parameters)
# solver.solve()
# in terminal do before u start:  
$ conda activate fenicsproject  
# Shortcuts:  
comment lines: cntr+K+C  
uncomment lines: cntr+K+U  
fold/unfold section: cntr+shift+[/]  
F8 Go to next error or warning  
Shift+F8 Go to previous error or warning  
go forward/back alt+arror  
  
# f_text = sympy.simplify(f_text) #make expression simpler then it is  
# x2 = 2*pi*pow(x, 3)+x*10+1/(x+1) #examples!  
# !!!calculate deriviation using sympy  
# Check additional solver parameters  
$ in console  
import fenics as f  
f.info(f.NonlinearVariationalSolver.default_parameters(), True)  
f.info(f.LinearVariationalSolver.default_parameters(), True)  
f.list_linear_solver_methods()  
  
У других солвером аналогично можно смотреть параметры!!!  
Также можно кастомизировать и другие параметры для конкретного солвера (там широкий спектр параметров для каждого)  
  
Jac     = derivative(T, w, dup) (пойдет?)  
# посмотреть содержимое пакета package content  
from fenics  import *
dir(package_name)  
$ for i in dir(numpy):  
$ print(i)  
# Junk from my file  
  
problem = NonlinearVariationalProblem(a == 0, u, bc)  
solver = NonlinearVariationalSolver(problem)  
solver_parameters={"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}  
solver.parameters.update(solver_parameters)  
solver.solve()  
# More junk  
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-10   # tolerance for coordinate comparisons
        return on_boundary and \
               (abs(x[1] - z1) < tol or abs(x[1] - z2) < tol or abs(x[0] - r2) < tol)s
# u0_boundary = DirichletBoundary()  
A1 = 2e-3  
A2 = 5e-3  
step = 1e-3  
array = numpy.linspace(A1, A2, 1+int((A2-A1)/step))  
  
if M is a class then M.__dict__  
# Move between cells  
cntrl+alt+[/]  
# Some junk
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
# Compare with analytical
fenics, numpy
vertex_values_u_D = interpolate(boundary.psi_sol_expr, V).compute_vertex_values(geometry.mesh)
vertex_values_u = u.compute_vertex_values(geometry.mesh)
error_max = numpy.max(numpy.abs(vertex_values_u_D - vertex_values_u))
print(error_max)

# Permeability class
class Permeability(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 1:
            values[0] = 1 # vessel
        else:
            values[0] = 0 # vacuum
    def value_shape(self):
        return ()
mu = Permeability(geometry.mesh, degree=0)

# get submesh of a mesh using index
submesh = SubMesh(geometry.mesh, subdomain_index)
matplt.xlim(geometry.r1, geometry.r2)
matplt.ylim(geometry.z1, geometry.z2)
plot(submesh)

# get u values on a submesh
submesh = SubMesh(geometry.mesh, subdomain_index)
u_plasma = u.compute_vertex_values(submesh)

# Newton solver part
a = dot(grad(u)/r, grad(r_2*v))*dx - etta * (p_coeff*r*r + F_2_coeff)*u*r*v*dx - tetta * sum(point_sources[2:len(point_sources)])*r*v*dx
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary)
solve(a == 0, u, bc)

# fenics version
$ dolfin-version 

# Expression with unknown function
funcc = Expression("u>=-0.02 ? 1 : 0", u=u, degree=2)
fu.countour_plot_via_mesh(geometry, interpolate(funcc, V), levels = p.levels, PATH = PATH, plot_title = '')

# Add my own module to anaconda path
conda develop ./Example/qwe.py -n fenicsproject
conda develop path/to/python/file -n fenicsproject
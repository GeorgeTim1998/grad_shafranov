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

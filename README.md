## RIGHT HAND SIDE f THAT YOU SEE IN PYTHON FILES IS INVERCED BECAUSE OF THE 
## VARIATIONSL FORMULATION  

## in terminal do before u start: 
# $source /home/george/anaconda3/bin/activate 
# $conda activate fenicsproject

## Shortcuts
# comment lines: cntr+K+C
# uncomment lines: cntr+K+U
# fold/unfold section: cntr+shift+[/]
# F8 Go to next error or warning
# Shift+F8 Go to previous error or warning
# go forward/back alt+arror 

{
    #r2weigth is r^2 that appears during
    #this function is used to define operator just like in G-Sh equation
    #make r2weigth available for spacial deriviations
    #if I want to state the problem as in
    # my "Вывод ур-ия Г-Ш.pdf" file then I probably need 'Expression' func
    # just like Ineed it to define right hand side of the equation
}

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

# problem = NonlinearVariationalProblem(a == 0, u, bc)
# solver = NonlinearVariationalSolver(problem)
# solver_parameters={"relative_tolerance": rel_tol, "absolute_tolerance": abs_tol}
# solver.parameters.update(solver_parameters)
# solver.solve()

More junk
# class DirichletBoundary(SubDomain):
#     def inside(self, x, on_boundary):
#         tol = 1E-10   # tolerance for coordinate comparisons
#         return on_boundary and \
#                (abs(x[1] - z1) < tol or abs(x[1] - z2) < tol or abs(x[0] - r2) < tol)

# u0_boundary = DirichletBoundary()

SQ_MIN = 2e-3
SQ_MAX = 5e-3
step = 1e-3
SQUARE_SIZE_ARRAY = numpy.linspace(SQ_MIN, SQ_MAX, 1+int((SQ_MAX-SQ_MIN)/step))
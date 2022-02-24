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

У других солвером аналогично можно смотреть параметры!!!
Также можно кастомизировать и другие параметры для конкретного солвера (там широкий спектр параметров для каждого)

Jac     = derivative(T, w, dup) (пойдет?)
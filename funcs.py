from __future__ import print_function
from fenics import *
import matplotlib.pyplot as matplt
# from matplotlib.pyplot import figure, gcf, isinteractive, title
import datetime
import numpy 
import sympy
from termcolor import colored
import pylab

DPI = 200

def Form_f_text(A1, A2):
    #A1 = 4*mo*p', A2 = FF'
    #deriviation are calculated using sympy library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("INVERCED right-hand equation side: \n", 'magenta') + f_text)

    return f_text

def Twod_plot(psi, x0, y1, y2, path): 
    # y1, y2 - min and max points in an interval of interest, 
    # x0 - point along which 2d graph is plotted
    # psi.set_allow_extrapolation(True)
    tol, point_num = 0.001, 100 + 1  # avoid hitting points outside the domain
    y = numpy.linspace(y1 + tol, y2 - tol, point_num)
    
    points = [(x0, y_) for y_ in y]  # create 2D points
    psi_line = numpy.array([psi(point) for point in points])
    
    matplt.plot(y, psi_line, 'k', linewidth=2)  # magnify w
    matplt.grid(True)
    matplt.xlabel('$r$')
    matplt.ylabel('$psi$')
    matplt.legend(["Point: %s, interval: [%s, %s]" % (x0, y1, y2), 'Load'], loc='upper left')
    
    time_title = Time_name()
    matplt.savefig('Figures/%s/%s_%s.png' % (path, time_title, x0), dpi = DPI)
    matplt.close() # close created plot
    
    print(colored("\t2d plot saved!", 'green'))
    
def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  #get current time to make figure name unique
    return time_title

def Save_figure(f_expr, mesh_r, mesh_z, addition, PATH):
    # move to funcs, add missing args, fix save path 
    # Plot solution and mesh. Save plot
    #nothing passed to function, because variables are global
    mesh_title = "%sx%s mesh" % (str(mesh_r), str(mesh_z))
    
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title) # file path+unique time name

    if addition == '_notitle':
        matplt.savefig("%s%s.png" %(path_my_file, addition), dpi = DPI) #no title figure for reports
    elif addition == '_title':
        matplt.title('Analyt Soloviev: %s\n%s' % (mesh_title, f_expr._cppcode)) # titled figure for my self
        matplt.savefig("%s%s.png" %(path_my_file, addition), dpi = DPI) #no title figure for reports
    else:
        matplt.title(addition) # titled figure for my self
        matplt.savefig("%s%s.png" %(path_my_file, addition), dpi = DPI) #no title figure for reports
    print(colored("\tPlot saved!", 'green'))
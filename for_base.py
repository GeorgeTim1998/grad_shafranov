import time
import numpy
import math
import datetime
import matplotlib.pyplot as matplt
from scipy.optimize import curve_fit

FONT = {'family': "Times New Roman",
        'size': 22}
matplt.rc('font', **FONT)
params = {'mathtext.default': 'regular' }          
matplt.rcParams.update(params)

def get_column(matrix, col):
    return [row[col] for row in matrix]

def func1(x, a, b):
    return a * x + b

def func2(x, a, b):
    return math.exp(b)*x**a


def plot_error_vs_mesh_from_file(folder_name, file_name, x_lim, PATH):
    with open("%s/%s.txt" % (folder_name, file_name), "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]

    mesh = get_column(data, 0)
    error = get_column(data, 1)
    
    popt, pcov = curve_fit(func1, numpy.log(mesh), numpy.log(error))    
    print(popt)
    matplt.loglog(numpy.array(mesh), func2(numpy.array(mesh), *popt), linewidth=1, label="$E_{max}$ ~ 1/$N_{\u2202\u03A9}^{2}$)")  # magnify w
    matplt.loglog(mesh, error, 'o', linewidth=2, label="$E_{max}$($N_{\u2202\u03A9}$)")  # magnify w
    matplt.legend()
    
    matplt.xlim(x_lim[0], x_lim[1])
    matplt.grid(True)
    matplt.xlabel("$N_{\u2202\u03A9}$")
    matplt.ylabel('$E_{max}$, Вб')

    save_contour_plot(PATH, '')

def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime) 
    return time_title

def save_contour_plot(PATH, plot_title):
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title)
    file_path = "%s.png" % path_my_file

    matplt.title(plot_title) 
    matplt.savefig(file_path, dpi=240, bbox_inches="tight")
    matplt.close()  

    print("3D countour plot saved to PATH: %s" % file_path)
    time.sleep(1)
    
plot_error_vs_mesh_from_file('Errors', 'Spheromak_09052022_224524', [0, 1050], "SHOW_NOW")
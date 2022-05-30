import time
import numpy
import math
import datetime
import matplotlib.pyplot as matplt
from scipy.optimize import curve_fit

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
    matplt.loglog(numpy.array(mesh), func2(numpy.array(mesh), *popt), linewidth=1)  # magnify w
    matplt.loglog(mesh, error, 'o', linewidth=2)  # magnify w
    
    matplt.xlim(x_lim[0], x_lim[1])
    matplt.grid(True)
    matplt.xlabel('Плотность сетки')
    matplt.ylabel('Максимальная ошибка, Вб')

    save_contour_plot(PATH, '')

def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  # get current time to make figure name unique
    return time_title

def save_contour_plot(PATH, plot_title):
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title)
    file_path = "%s.png" % path_my_file

    matplt.title(plot_title)  # titled figure for my self
    # no title figure for reports
    matplt.savefig(file_path, dpi=240, bbox_inches="tight")
    matplt.close()  # close created plot

    print("3D countour plot saved to PATH: %s" % file_path)
    time.sleep(1)
    
plot_error_vs_mesh_from_file('Errors', '1D_problem_09052022_235717', [0, 1000], "SHOW_NOW")
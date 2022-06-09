# import mshr
import fenics as fn

import numpy as np
import matplotlib.pyplot as matplt
import matplotlib.tri as tri
import funcs as fu
from geometry import Geometry
import point_source_data as psd
import numpy as np

# ge = Geometry()
# ps = psd.PointSource(1)
# r = np.transpose(ps.r)[0]
# z = np.transpose(ps.r)[1]

# ge.outer_mephist_vessel()
# ge.inner_mephist_vessel()
# matplt.plot(ge.outer_vessel_contour[0],
#             ge.outer_vessel_contour[1],
#             c='k', linewidth=1)
# matplt.plot(ge.inner_vessel_contour[0],
#             ge.inner_vessel_contour[1],
#             c='k', linewidth=1)

# matplt.scatter(r, z, c='g')

# matplt.xlabel("r, м")
# matplt.ylabel("z, м")
# matplt.xticks([0.1, 0.2, 0.3, 0.4, 0.5])

# matplt.xlim([0.05, 0.55])
# matplt.ylim([-0.4, 0.4])
# # matplt.grid(True)
# matplt.gca().set_aspect("equal")
# fu.save_contour_plot(PATH="SHOW_NOW", plot_title='')

fu.plot_error_vs_mesh_from_file('Errors', 'Spheromak_09052022_224524', [0, 1000], "SHOW_NOW")

# fu.shrink_contour('Data',
#                   'Mephist_vessel_outer_surface',
#                   shrunk_to_point = [0.25, 0],
#                   alpha_x=0.7,
#                   alpha_z=0.8)

# def get_column(matrix, col):
#     return [row[col] for row in matrix]
    
# def read_data_from_file():
#     file_path = "Data/Mephist_vessel_outer_surface.txt"
#     with open(file_path, "r") as file:  # change to Read_from_file func
#         data = [[float(num) for num in line.split(',')] for line in file]

#     x = np.array(get_column(data, 0))
#     z = np.array(get_column(data, 1))

#     return x, z

# [x,z] = read_data_from_file()
    
# matplt.plot(x, z)
# matplt.xticks(np.linspace(0.1, 0.4, 7))
# matplt.grid(True)
# matplt.gca().set_aspect("equal")
# # matplt.savefig('Data/Mephist_vessel_outer_surface.svg',dpi=350)
# matplt.show()

# mesh = fn.RectangleMesh(fn.Point(-1, -1), fn.Point(1, 1), 100, 100)
# f = fn.Expression("10*(1-pow(x[0] - 1, 2)-pow(x[1]-1, 2))", degree=2)
# V = fn.FunctionSpace(mesh, 'P', 1)
# print(fn.interpolate(f, V).vector().max())
# fn.plot(fn.interpolate(f, V))
# matplt.show()

# fu.plot_error_vs_mesh_from_file(folder_name='Errors', file_name='1D_problem_09052022_235717', x_lim=[0, 1000], PATH='1D_problem')
# fu.plot_Dina_results("DINA")

# x, y, z = np.genfromtxt(r'psi.dat', unpack=True)

# # x = np.unique(x)
# # y = np.unique(y)
# # z = z.reshape((len(y), len(x)))

# # triang = tri.Triangulation(x, y)
# # interpolator = tri.LinearTriInterpolator(triang, z)

# # Xi, Yi = np.meshgrid(x, y)
# # zi = interpolator(Xi, Yi)

# matplt.tricontour(x, y, z, levels = 100)
# matplt.colorbar()
# matplt.gca().set_aspect("equal")
# matplt.show()

# def Column(matrix, col):
#     return [row[col] for row in matrix]

# file_path = "КонтурКамерыСредний.txt"
# with open(file_path, "r") as file: # change to Read_from_file func
#     data = [[float(num) for num in line.split('  ')] for line in file]

# x = numpy.array(Column(data, 0))*1e-3
# z = numpy.array(Column(data, 1))*1e-3

# camera_height = 581.69*1e-3
# z_height = z[-2]

# x = x[0:len(x)-1]
# z = z[0:len(z)-1]

# z = z - z_height/2 * numpy.ones(len(z))
# z = z + (camera_height/2 - numpy.amax(z)) * numpy.ones(len(z)) 
# x = numpy.append(x, numpy.flip(x)) # move it along x axis!
# z = numpy.append(z, numpy.flip(-z))

# print(numpy.amax(x))
# print(numpy.amin(x))
# print(numpy.amax(z))

# print(numpy.amin(z))

# matplt.plot(x, z)
# matplt.grid(true)
# matplt.gca().set_aspect("equal")
# matplt.xlim(0, 0.6)
# matplt.ylim(-0.4, 0.4)
# matplt.show()

# point_list = []

# for i in range(len(x)):
#     point_list.append(fenics.Point(x[i], z[i]))
# # a = fenics.Point(0,0)
# # b = fenics.Point(4,0)
# # c = fenics.Point(10,0.7)
# # d = fenics.Point(2,5)

# domain = mshr.Polygon(point_list)
# mesh = mshr.generate_mesh(domain, 20)

# fenics.plot(mesh)
# matplt.show()

# # *** Error:   Unable to create polygon.
# # *** Reason:  Polygon vertices must be given in counter clockwise order.
# # *** Where:   This error was encountered inside CSGPrimitives2D.cpp.
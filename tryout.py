# import mshr
# import fenics

import numpy as np
import matplotlib.pyplot as matplt
import matplotlib.tri as tri
import funcs as fu

fu.plot_Dina_results("DINA")

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
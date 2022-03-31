import mshr
import fenics
import matplotlib.pyplot as matplt

a = fenics.Point(0,0)
b = fenics.Point(4,0)
c = fenics.Point(10,0.7)
d = fenics.Point(2,5)

domain = mshr.Polygon([a,b,c,d])
mesh = mshr.generate_mesh(domain, 20)

fenics.plot(mesh)
matplt.show()

# *** Error:   Unable to create polygon.
# *** Reason:  Polygon vertices must be given in counter clockwise order.
# *** Where:   This error was encountered inside CSGPrimitives2D.cpp.
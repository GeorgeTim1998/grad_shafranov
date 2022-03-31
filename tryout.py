import mshr
import fenics
import matplotlib.pyplot as matplt

a = fenics.Point(0,0)
b = fenics.Point(1,0)
c = fenics.Point(0.5,0.5)

domain = mshr.Polygon([a,b,c])
mesh = mshr.generate_mesh(domain, 10)

fenics.plot(domain)
matplt.show()

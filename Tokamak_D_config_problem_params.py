import numpy
import logger

class Problem:
    def __init__(self):
        self.mesh_density = 50
        self.default_mesh = 100
        self.default_mesh_array = numpy.linspace(50, 1000, 1+int((1000-50)/50))
        self.domain_boundary_coordinates = [0.02, 1.5, 
                                            -0.7, 0.7]
        
        self.psi_0 = 0.1
        self.eps = 0.05
        self.a = 1
        self.b = 0.5
        
        self.levels = 20
        self.levels = numpy.linspace(-0.10, 10, 25)
        
        self.errors = []
        
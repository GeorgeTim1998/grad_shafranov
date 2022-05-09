import numpy
import logger

class Problem:
    def __init__(self):
        self.mesh_density = 50
        self.default_mesh = 80
        self.default_mesh_array = numpy.linspace(50, 1000, 1+int((1000-50)/50))
        self.domain_boundary_coordinates = [0.05, 0.7, 
                                            -0.6, 0.6]
        
        self.psi_0 = 0.1
        self.R = 0.3
        self.alpha = 1.1
        
        self.levels = 20
        self.levels = numpy.linspace(-0.14, self.psi_0, 25)
        
        self.errors = []
        
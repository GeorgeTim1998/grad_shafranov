import numpy
import logger

class Problem:
    def __init__(self):
        self.mesh_density = 50
        self.default_mesh = 500
        self.default_mesh_array = numpy.linspace(50, 950, 1+int((950-50)/50))
        self.domain_boundary_coordinates = [0, 1.25, 
                                            -0.7, 0.7]
        
        self.psi_0 = 0.1
        self.eps = 0.05
        self.a = 1
        self.b = 0.4
        
        self.psi_min = -1
        self.psi_max = 9
        self.psi_step = 1
        
        self.levels = 20
        self.levels = numpy.linspace(self.psi_min, self.psi_max, 1+int((self.psi_max-self.psi_min)/self.psi_step))
        
        self.errors = []
        
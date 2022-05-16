from cmath import sqrt
import numpy
import logger

class Problem:
    def __init__(self):
        self.mesh_density = 50
        self.default_mesh = 150
        self.default_mesh_array = numpy.linspace(50, 950, 1+int((950-50)/50))
        self.domain_boundary_coordinates = [0, 0.5, 
                                            -0.6, 0.6]
        
        self.psi_0 = 0.1
        self.R = 0.2*sqrt(2)
        self.alpha = 0.5
        
        self.psi_min = -0.02
        self.psi_max = self.psi_0
        self.psi_step = 0.01
        
        self.levels = 20
        self.levels = numpy.linspace(self.psi_min, self.psi_max, 1+int((self.psi_max-self.psi_min)/self.psi_step))
        
        self.errors = []
        
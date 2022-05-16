import numpy
import logger

class Problem:
    def __init__(self):
        self.default_mesh = 500
        self.default_mesh_array = numpy.linspace(50, 950, 1+int((950-50)/50))
        self.domain_boundary_coordinates = [0.02, 1.5, 
                                            -0.7, 0.7]
        
        self.a = 0.1
        self.b = 2
        self.j0 = 1
        
        self.errors = []
        
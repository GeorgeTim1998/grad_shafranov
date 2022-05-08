#%% Imports
from fenics import *
import matplotlib.pyplot as matplt
import logger
import mshr
import time
import funcs as fu
import MEPHIST_data as M
import logger
from geometry import Geometry
from boundary_conditions import BoundaryConditions
import numpy

class Problem:
    def __init__(self):
        self.domain_geometry = [0.01, 1.04, -0.5, 0.5]
        
        self.default_mesh = 500
        
        self.boundary_condition_str = '-1e-5'
        
        self.m0_p2 = 0.277
        self.psi_0 = 0.1
        
        self.A1 = 1e-3  
        self.A2 = 1e-2  
        self.step = 0.5e-3
        self.cycle_array = numpy.linspace(self.A1, self.A2, 1+int((self.A2-self.A1)/self.step))  
            
        self.contour_levels = 20
        # self.contour_levels = list(numpy.geomspace(-1e-3, 1e-8))
import numpy
import MEPHIST_data as M
import logger
class Problem:
    def __init__(self):
        self.domain_geometry = [0.05, 0.55, -0.4, 0.4]
        
        self.plasma_centre_point = [0.235, 0]
        self.plasma_radius = 0.09
        self.plasma_domain_segments=100
        self.point_source_disp = 1
        
        self.mesh_density = 180
        
        self.boundary_condition_str = '0'
        
        self.betta = 0.1
        self.psi_axis = 0.05777839
        self.psi_pl_edge = 0.015
        self.psi_pl_expr = 0.007
        
        self.p_correction = 1
        self.F_correction = 1
        self.psi_correction = 1
        
        self.p_correction_array = numpy.geomspace(0.1, 10, 5)
        self.F_correction_array = numpy.geomspace(0.1, 10, 5)
        self.psi_correction_array = numpy.geomspace(0.1, 10, 5)
        
        self.A1 = 1e-3  
        self.A2 = 1e-2  
        self.step = 0.5e-3
        self.cycle_array = numpy.linspace(self.A1, self.A2, 1+int((self.A2-self.A1)/self.step))  
            
        self.levels = 10
        # self.levels = numpy.array([0.007,0.0072, 0.0074, 0.0076, 0.0078, 0.008])
        self.levels = numpy.array([-0.0007, -0.0004, -0.0002, 0, 0.001, 0.003,0.005, 0.007, 0.01])
        # self.levels = numpy.linspace(-0.0007, 0.01, )
        logger.log_n_output_colored_message(colored_message="point_source_disp = ", color='green', white_message=str(self.point_source_disp))
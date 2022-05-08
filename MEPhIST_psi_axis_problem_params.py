import numpy
import MEPHIST_data as M
import logger
class Problem:
    def __init__(self):
        self.domain_geometry = [0.05, 0.55, -0.4, 0.4]
        
        self.plasma_centre_point = [0.23, 0]
        self.plasma_radius = M.MEPhIST().a*0.65
        self.plasma_domain_segments=60
        self.point_source_disp = 1
        
        self.mesh_density = 180
        
        self.boundary_condition_str = '0'
        
        self.p_correction = 1
        self.F_correction = 1
        self.psi_correction = 1
        
        self.A1 = 1e-3  
        self.A2 = 1e-2  
        self.step = 0.5e-3
        self.cycle_array = numpy.linspace(self.A1, self.A2, 1+int((self.A2-self.A1)/self.step))  
            
        self.contour_levels = 150
        # self.contour_levels = list(numpy.geomspace(-1e-3, 1e-8))
        logger.log_n_output_colored_message(colored_message="point_source_disp = ", color='green', white_message=str(self.point_source_disp))
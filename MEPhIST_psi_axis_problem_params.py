import MEPHIST_data as M
class Problem:
    def __init__(self):
        self.domain_geometry = [0.05, 0.55, -0.4, 0.4]
        
        self.plasma_centre_point = [0.23, 0]
        self.plasma_radius = M.MEPhIST().a*0.65
        self.plasma_domain_segments=60
        
        self.mesh_density = 180
        
        self.boundary_condition_str = '-1e-2'
        
        self.p_correction = 1
        self.F_correction = 1
        self.psi_correction = 1
        
        self.contour_levels = 40
        # levels = list(numpy.geomspace(-1e-3, 1e-8))
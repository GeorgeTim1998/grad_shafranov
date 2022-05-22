import numpy
import MEPHIST_data as M
import logger
import math
import logger


class Problem:
    def __init__(self):
        self.domain_geometry = [0.05, 0.55, -0.4, 0.4]
        self.mesh_density = 120

#%% Vessel info
        self.centre_point = [0.3, 0]
        self.ves_inner_radius = 0.1  # inner boundary info
        self.ves_inner_segments = 60
        
        self.ves_outer_radius = 0.15  # outer boundary info
        self.ves_outer_segments = 100
        
        self.surface_ratio = math.pi
#%% Physical properties
        self.VACUUM_PERMEABILITY = 1
        self.VESSEL_PERMEABILITY = 1.008
        self.PLASMA_PERMEABILITY = 1

        self.VACUUM_CONDUCTIVITY = 0
        self.VESSEL_CONDUCTIVITY = 7.2e-6
        self.PLASMA_CONDUCTIVITY = 1e-6
#%% Problem params
        self.psi_0 = 1e-3
        self.R = 0.2*math.sqrt(2)
        self.alpha = 0.5
        
        self.boundary_condition_str = '0'
#%% Plotting arrays and levels
        self.p_correction_array = numpy.geomspace(0.1, 10, 5)
        self.F_correction_array = numpy.geomspace(0.1, 10, 5)
        self.psi_correction_array = numpy.geomspace(0.1, 10, 5)

        self.A1 = 1e-3
        self.A2 = 1e-2
        self.step = 0.5e-3
        self.cycle_array = numpy.linspace(
            self.A1, self.A2, 1+int((self.A2-self.A1)/self.step))

        self.levels = 20
        # self.levels = numpy.array([0.007,0.0072, 0.0074, 0.0076, 0.0078, 0.008])
        # self.levels = numpy.array([-0.0007, -0.0004, -0.0002, 0, 0.001, 0.003,0.005, 0.007, 0.01])
        # self.levels = numpy.linspace(-0.0007, 0.01, )

        logger.log_n_output_colored_message(
            colored_message="VACUUM_PERMEABILITY = ", color='green', white_message=str(self.VACUUM_PERMEABILITY))
        logger.log_n_output_colored_message(
            colored_message="VESSEL_PERMEABILITY = ", color='green', white_message=str(self.VESSEL_PERMEABILITY))
        
        logger.log_n_output_colored_message(
            colored_message="VESSEL_CONDUCTIVITY = ", color='green', white_message=str(self.VESSEL_CONDUCTIVITY))
        logger.log_n_output_colored_message(
            colored_message="PLASMA_CONDUCTIVITY = ", color='green', white_message=str(self.PLASMA_CONDUCTIVITY))

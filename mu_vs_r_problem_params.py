import numpy
import MEPHIST_data as M
import logger
import math
import logger
from funcs import M0


class Problem:
    def __init__(self):
        self.domain_geometry = [0.05, 0.55, -0.4, 0.4]
        self.mesh_density = 120

# %% Vessel info
        self.centre_point = [0.3, 0]
        self.ves_inner_radius = 0.1  # inner boundary info
        self.ves_inner_segments = 100

        self.ves_outer_radius = 0.15  # outer boundary info
        self.ves_outer_segments = int(
            self.ves_outer_radius/self.ves_inner_radius * self.ves_inner_segments)

        self.area_ratio = 100*math.pi*(self.ves_outer_radius**2 - self.ves_inner_radius**2) / (
            self.domain_geometry[1]-self.domain_geometry[0]) / (self.domain_geometry[3]-self.domain_geometry[2])
# %% Physical properties
        self.VACUUM_PERMEABILITY = 1
        self.VESSEL_PERMEABILITY = 1.008
        self.PLASMA_PERMEABILITY = 1

        self.VACUUM_CONDUCTIVITY = 0
        self.VESSEL_CONDUCTIVITY = 1.4e6
        self.PLASMA_CONDUCTIVITY = 9e6
# %% Problem params
        self.psi_0 = 1e-3
        self.R = 0.3
        self.alpha = 0.5

        self.ts = M0 * self.VESSEL_PERMEABILITY * self.VESSEL_CONDUCTIVITY * \
            (self.ves_outer_radius - self.ves_inner_radius)**2
            
        self.t0 = 0
        self.t_max = 1e-1*self.ts
        num_of_t = 2+1
        self.t = numpy.linspace(self.t0, self.t_max, num_of_t)

        self.tm = 0.5*1e-1*self.ts
        
        self.disp_fact = 1 # множитель для характерного смещения на радиус камеры

        self.boundary_condition_str = '0'
# %% Plotting arrays and levels
        self.p_correction_array = numpy.geomspace(0.1, 10, 5)
        self.F_correction_array = numpy.geomspace(0.1, 10, 5)
        self.psi_correction_array = numpy.geomspace(0.1, 10, 5)

        self.A1 = 1e-3
        self.A2 = 1e-2
        self.step = 0.5e-3
        self.cycle_array = numpy.linspace(
            self.A1, self.A2, 1+int((self.A2-self.A1)/self.step))

        self.levels = 20
        self.amount_of_levels = 20 # for a method below! (find_levels)
        # self.levels = numpy.array([0.007,0.0072, 0.0074, 0.0076, 0.0078, 0.008])
        # self.levels = numpy.array([-0.0007, -0.0004, -0.0002, 0, 0.001, 0.003,0.005, 0.007, 0.01])
        # self.levels = numpy.linspace(-0.0007, 0.01, )
# %% Logs
        logger.info(message="t = linspace(%.4e, %.4e, %d)" %
                    (self.t0, self.t_max, int(num_of_t)))
        logger.log_n_output_colored_message(
            colored_message="ts = ", color='green', white_message=str(self.ts)
        )
        logger.log_n_output_colored_message(
            colored_message="tm = ", color='green', white_message=str(self.tm)
        )

        logger.log_n_output_colored_message(
            colored_message="VACUUM_PERMEABILITY = ", color='green', white_message=str(self.VACUUM_PERMEABILITY))
        logger.log_n_output_colored_message(
            colored_message="VESSEL_PERMEABILITY = ", color='green', white_message=str(self.VESSEL_PERMEABILITY))
        logger.log_n_output_colored_message(
            colored_message="PLASMA_PERMEABILITY = ", color='green', white_message=str(self.PLASMA_PERMEABILITY))

        logger.log_n_output_colored_message(
            colored_message="VACUUM_CONDUCTIVITY = ", color='green', white_message=str(self.VACUUM_CONDUCTIVITY))
        logger.log_n_output_colored_message(
            colored_message="VESSEL_CONDUCTIVITY = ", color='green', white_message=str(self.VESSEL_CONDUCTIVITY))
        logger.log_n_output_colored_message(
            colored_message="PLASMA_CONDUCTIVITY = ", color='green', white_message=str(self.PLASMA_CONDUCTIVITY))

        logger.log_n_output_colored_message(
            colored_message="Area ratio: ", color='green', white_message=str(self.area_ratio)
        )

#%% Find levels
    def find_levels(self, u):
        u_min = 0
        u_max = u_max = u.vector()[:].max()
        
        exponent = math.ceil(abs(math.log10(u_max))) + 1
        u_max = int(u_max * pow(10, exponent)) / pow(10, exponent)
        self.levels = numpy.linspace(u_min, u_max, 11)
    

import numpy
import MEPHIST_data as M
import logger
import math
import logger
from funcs import M0


class Problem:
    def __init__(self):
        self.domain_factors = numpy.array([1, 3, 2.5, 2.5])#numpy.array([1, 2, 1.5, 1.5])
        self.domain_geometry0 = numpy.array([0.05, 13, -9, 9])
        self.domain_geometry = self.domain_factors * self.domain_geometry0
        
        self.plot_domain = self.domain_geometry0
        self.mesh_density = 240

        # Vessel info
        self.centre_point = [0.2, 0]
        
        self.vessel_thickness = 0.05
        self.plasma_step_length = 0.02

        # Physical properties of plasma and vacuum vessel
        self.VACUUM_PERMEABILITY = 1
        self.VESSEL_PERMEABILITY = 1.008
        self.PLASMA_PERMEABILITY = 1

        self.VACUUM_CONDUCTIVITY = 0
        self.VESSEL_CONDUCTIVITY = 1.4e7 #1.4e6
        self.PLASMA_CONDUCTIVITY = 0 #8e5 #1.4e5 # 9e6 - посчитанно с использованием montani2021: sigma ~ T**3/2

        # Problem params
        self.R = self.centre_point[0]

        self.ts = M0 \
            * self.VESSEL_PERMEABILITY * self.VESSEL_CONDUCTIVITY \
            * (self.vessel_thickness)**2

        self.t0 = 0
        self.ts_fraction = 0.5
        self.ts_fraction_tm = 0.1
        self.t_max = self.ts_fraction *0.004410000000000001#self.ts_fraction * self.ts

        self.disp_fact = 1  # множитель характерного смещения. Умножаем на радиус камеры
        self.tm = self.ts_fraction_tm * 0.004410000000000001#self.ts

        self.num_of_t = 1  # число точек по времени +1 из-за 0
        self.t = numpy.linspace(self.t0, self.t_max,
                                self.num_of_t)  # Временной массив

        # Характерное время смещения. Вычисляется как часть скинового времени

        self.boundary_condition_str = '0'

        # Plotting arrays and levels
        self.levels = 50
        self.amount_of_levels = 20  # for a method below! (find_levels)
        self.step = 0.1  # Шаг для линий уровня. Если максимум ф-ии 1.2e-3, то шаг будет self.step*e-3

        # self.__log_problem_params()

    def find_levels(self, u, step):  # Find levels for plotting
        u_min = u.vector()[:].min()
        u_max = u.vector()[:].max()

        exponent = math.ceil(
            abs(math.log10(abs(u_max))))# u_max
        # Найти ближайшее снизу к u_max число, чтобы отличие было максимум на величину шага

        u_max = step \
            * pow(10, -exponent) \
            * int((u_max) / (step*pow(10, -exponent)))
        u_min = step \
            * pow(10, -exponent) \
            * int((u_min) / (step*pow(10, -exponent)))

        self.levels = numpy.linspace(u_min, u_max, int(
            (u_max-u_min) / (step*pow(10, -exponent))) + 1
        )
        
    def find_levels_with_exponent(self, u, exponent, step):
        u_min = 0
        u_max = u.vector()[:].max()

        u_exp = math.ceil(abs(math.log10(u_max)))
        u_max = int(
            u_max * pow(10, u_exp)
        ) / pow(10, u_exp)

        self.levels = numpy.linspace(u_min, u_max, int(
            (u_max - u_min)/(step * pow(10, -exponent))) + 1
        )
        
    # def __log_problem_params(self):
    #     logger.info(message="t = linspace(%.4e, %.4e, %d)" %
    #                 (self.t0, self.t_max, int(self.num_of_t)))
    #     logger.log_n_output_colored_message(
    #         colored_message="ts = ",
    #         color='green',
    #         white_message=str(self.ts)
    #     )
    #     logger.log_n_output_colored_message(
    #         colored_message="tm = ",
    #         color='green',
    #         white_message=str(self.tm)
    #     )
    #     logger.log_n_output_colored_message(
    #         colored_message="ts_fract_t_max = ",
    #         color='green',
    #         white_message=str(self.ts_fraction)
    #     )
    #     logger.log_n_output_colored_message(
    #         colored_message="ts_fraction_tm = ",
    #         color='green',
    #         white_message=str(self.ts_fraction_tm)
    #     )

    #     logger.log_n_output_colored_message(
    #         colored_message="VACUUM_PERMEABILITY = ", 
    #         color='green', 
    #         white_message=str(self.VACUUM_PERMEABILITY))
    #     logger.log_n_output_colored_message(
    #         colored_message="VESSEL_PERMEABILITY = ", 
    #         color='green', 
    #         white_message=str(self.VESSEL_PERMEABILITY))
    #     logger.log_n_output_colored_message(
    #         colored_message="PLASMA_PERMEABILITY = ", 
    #         color='green', 
    #         white_message=str(self.PLASMA_PERMEABILITY))

    #     logger.log_n_output_colored_message(
    #         colored_message="VACUUM_CONDUCTIVITY = ", 
    #         color='green', 
    #         white_message=str(self.VACUUM_CONDUCTIVITY))
    #     logger.log_n_output_colored_message(
    #         colored_message="VESSEL_CONDUCTIVITY = ", 
    #         color='green', 
    #         white_message=str(self.VESSEL_CONDUCTIVITY))
    #     logger.log_n_output_colored_message(
    #         colored_message="PLASMA_CONDUCTIVITY = ", 
    #         color='green', 
    #         white_message=str(self.PLASMA_CONDUCTIVITY))



    def __find_exponent(self, u):
        u_max = u.vector()[:].max()
        self.exponent = math.ceil(abs(math.log10(u_max)))

    # def __calucate_domain_area(self):
    #     delta_r = self.domain_geometry[1]-self.domain_geometry[0]
    #     delta_z = self.domain_geometry[3]-self.domain_geometry[2]

    #     return delta_r*delta_z

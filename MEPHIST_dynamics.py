# Решаем динамическую задачу движения точечного источника в камере МИФИСТа, вместе со всей полоидальной системой МИФИСТа

# %% Imports
import sys
import time
import numpy
import matplotlib.pyplot as matplt

# from . boundary_conditions import BoundaryConditions
# from .. import boundary_conditions
# from expressions import Expressions
# from geometry import Geometry
# # import MEPHIST_dynamics_params as Problem
# import point_source_data as psd
# import funcs as fu
# import logger

from fenics import *


# %% Pre-programm stuff
PATH = 'MEPHIST_dynamics'

t0 = time.time()
# current_pyfile = "\n\n---------MEPHIST_dynamics.py---------"
# logger.log_n_output("%s" % current_pyfile, 'red')

# # %% Needed objects and contour levels
# boundary_conditions = boundary_conditions.BoundaryConditions()
# geometry = Geometry()
# # p = Problem()
# e = Expressions()

# # %% Domain and mesh definition
# # domain = geometry.rectangle_domain(
# #     area=[p.domain_geometry[0], p.domain_geometry[1], p.domain_geometry[2], p.domain_geometry[3]])

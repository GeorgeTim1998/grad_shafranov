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
import point_source_data as psd
import numpy
import MEPhIST_2_problems_problem_params as P
    
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------MEPhIST_2_problems.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'MEPhIST_2_problems'

#%% Needed objects
boundary_conditions = BoundaryConditions()
geometry = Geometry()
psi_axis = M.MEPhIST().psi_axis
problem = P.Problem()

#%% Domain and mesh definition
domain = geometry.circle_domain(centre_point=problem.plasma_centre_point, radius=problem.plasma_radius, segments=problem.plasma_domain_segments)
geometry.generate_mesh_in_domain(domain=domain, density=problem.mesh_density)

#%% Define function space
V = FunctionSpace(geometry.mesh, 'Lagrange', 1)

u = TrialFunction(V)
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition(problem.boundary_condition_str)
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary)

#%% Solve
[r_2, r] = geometry.operator_weights(V)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(problem.point_source_disp)))

for psi_correction in problem.psi_correction_array:
    for F_correction in problem.F_correction_array:
        for p_correction in problem.p_correction_array:
            [p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_pow_2_iteration(p_correction=p_correction, F_correction=F_correction, psi_axis=psi_correction*problem.psi_axis)
            logger.log_n_output_colored_message(colored_message="Correction coeff for psi on axis = ", color='green', white_message=str(problem.psi_correction))

            u = TrialFunction(V)
            a = dot(grad(u)/r, grad(r_2*v))*dx - (p_coeff*r*r + F_2_coeff)*u*r*v*dx
            L = Constant(0)*r*v*dx
            # L = tetta * sum(point_sources[2:len(point_sources)])*r*v*dx

            u = Function(V)
            solve(a == L, u, bc)

            #%% Post solve
            fu.What_time_is_it(t0, 'Variational problem solved')
            fu.countour_plot_via_mesh(geometry, u, levels = problem.contour_levels, PATH = PATH, plot_title = '')

            # fu.fenics_plot(u, PATH, plot_title='')

            fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
            logger.log_n_output_colored_message(colored_message="'Done'\n", color='red', white_message='')
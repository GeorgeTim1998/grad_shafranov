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
import MEPhIST_psi_axis_problem_params as P
    
#%% Pre-programm stuff
t0 = time.time()
current_pyfile = '---------MEPhIST_psi_axis.py---------'
logger.log_n_output("%s" % current_pyfile, 'red')
fu.print_colored("Date_Time is: %s" % fu.Time_name(), 'cyan')
PATH = 'MEPhIST_psi_axis'

#%% Needed objects and contour levels
boundary_conditions = BoundaryConditions()
geometry = Geometry()
psi_axis = M.MEPhIST().psi_axis
problem = P.Problem()

#%% Domain and mesh definition
domain = geometry.rectangle_domain(area=problem.domain_geometry)
plasma_circle = geometry.circle_domain(centre_point=problem.plasma_centre_point, radius=problem.plasma_radius, segments=problem.plasma_domain_segments)

domain.set_subdomain(1, plasma_circle)

geometry.generate_mesh_in_domain(domain=domain, density=problem.mesh_density)

markers = MeshFunction("size_t", geometry.mesh, geometry.mesh.topology().dim(), geometry.mesh.domains())

#%% Define function space and step coefficients
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh

u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition(problem.boundary_condition_str)
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(1)))

# A1 = 1e-3  
# A2 = 1e-2  
# step = 0.5e-3
# array = numpy.linspace(A1, A2, 1+int((A2-A1)/step))  

# for correction in array:
[p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_pow_2_iteration(p_correction=problem.p_correction, F_correction=problem.F_correction, psi_axis=problem.psi_correction*psi_axis)
logger.log_n_output_colored_message(colored_message="Correction coeff for psi on axis = ", color='green', white_message=str(problem.psi_correction))

u = TrialFunction(V)
a = dot(grad(u)/r, grad(r_2*v))*dx - (p_coeff*r*r + F_2_coeff)*u*r*v*dx(1)
L = sum(point_sources)*r*v*dx(0)

u = Function(V)
solve(a == L, u, bc)

#%% Post solve
fu.What_time_is_it(t0, 'Variational problem solved')
fu.countour_plot_via_mesh(geometry, u, levels = problem.contour_levels, PATH = PATH, plot_title = '')

fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
logger.log_n_output_colored_message(colored_message="'Done'\n", color='red', white_message='')
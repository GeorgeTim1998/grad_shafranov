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
import MEPHIST_data as MEPH
    
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

levels = 40
# levels = list(numpy.geomspace(-1e-3, 1e-8))

#%% Domain and mesh definition
domain = geometry.rectangle_domain(area=[0.05, 0.55, -0.4, 0.4])
mephist_vessel = geometry.mephist_vessel()
plasma_circle = geometry.circle_domain(centre_point=[0.23, 0], radius=M.MEPhIST().a*0.65, segments=60)

no_plasma_domain = mephist_vessel - plasma_circle

domain.set_subdomain(1, no_plasma_domain)
domain.set_subdomain(2, plasma_circle)

geometry.generate_mesh_in_domain(domain=domain, density=100)

markers = MeshFunction("size_t", geometry.mesh, geometry.mesh.topology().dim(), geometry.mesh.domains())

#%% Define function space and step coefficients
V = FunctionSpace(geometry.mesh, 'P', 1) # standard triangular mesh

u = TrialFunction(V) # u must be defined as function before expression def
v = TestFunction(V)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition("0")
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary) #гран условие как в задаче дирихле

#%% Solve
[r_2, r] = geometry.operator_weights(V)

dx = Measure('dx', domain=geometry.mesh, subdomain_data=markers)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(1)))

A1 = 1e-3  
A2 = 1e-2  
step = 0.5e-3
array = numpy.linspace(A1, A2, 1+int((A2-A1)/step))  

for correction in array:
    [p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_pow_2_iteration(p_correction=1, F_correction=1, psi_axis=correction*psi_axis)
    logger.log_n_output_colored_message(colored_message="correction=", color='red', white_message=str(correction))
    
    u = TrialFunction(V)
    a = dot(grad(u)/r, grad(r_2*v))*dx - (p_coeff*r*r + F_2_coeff)*u*r*v*dx(2)
    L = (sum(point_sources)*r*v*dx(0) + sum(point_sources)*r*v*dx(1))

    u = Function(V)
    solve(a == L, u, bc)

    #%% Post solve
    fu.What_time_is_it(t0, 'Variational problem solved')
    fu.countour_plot_via_mesh(geometry, u, levels = levels, PATH = PATH, plot_title = '')

    fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
    logger.log_n_output_colored_message(colored_message="'Done'\n", color='red', white_message='')
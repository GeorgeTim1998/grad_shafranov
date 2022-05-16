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
import math as m
    
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
p = P.Problem()

#%% Domain and mesh definition
domain = geometry.rectangle_domain(area=p.domain_geometry)
plasma_circle = geometry.circle_domain(centre_point=p.plasma_centre_point, radius=p.plasma_radius, segments=p.plasma_domain_segments)
domain.set_subdomain(1, plasma_circle)
# domain = geometry.circle_domain(centre_point=p.plasma_centre_point, radius=p.plasma_radius, segments=p.plasma_domain_segments)
geometry.generate_mesh_in_domain(domain=domain, density=p.mesh_density)

markers = MeshFunction("size_t", geometry.mesh, geometry.mesh.topology().dim(), geometry.mesh.domains())

class PlasmaStepConstant(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 1:
            values[0] = 1 # vessel
        else:
            values[0] = 0 # vacuum
    def value_shape(self):
        return ()

class VacuumStepConstant(UserExpression):
    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 1:
            values[0] = 0 # vessel
        else:
            values[0] = 1 # vacuum
    def value_shape(self):
        return ()
    
etta = PlasmaStepConstant(geometry.mesh, degree=1)
tetta = VacuumStepConstant(geometry.mesh, degree=1)
#%% Define function space
V = FunctionSpace(geometry.mesh, 'Lagrange', 1)

#%% Boundary conditions
u_D = boundary_conditions.constant_boundary_condition(p.boundary_condition_str)
bc = DirichletBC(V, u_D, fu.Dirichlet_boundary)

#%% Solve
[r_2, r] = geometry.operator_weights(V)

point_sources = fu.Array_Expression(fu.ArrayOfPointSources(psd.PointSource(p.point_source_disp)))

[p_coeff, F_2_coeff] = fu.plasma_sources_coefficients_exp_profile(p_correction=p.p_correction, F_correction=p.F_correction, psi_correction=p.psi_correction)

u = Function(V)
v = TestFunction(V)

po_2 = (u-Constant(p.psi_pl_edge))/Constant(p.psi_axis-p.psi_pl_edge)
G_psi = (exp(1 - po_2) - 1)/(m.e - 1)

L = tetta * sum(point_sources[2:len(point_sources)])*r*v*dx
a = dot(grad(u)/r, grad(r_2*v))*dx + etta*(p_coeff*r*r + F_2_coeff) * G_psi * r*v*dx - L

du = TrialFunction(V)
J = derivative(a, u, du)
solve(a == 0, u, bc, J=J)

#%% Post solve
fu.What_time_is_it(t0, 'Variational problem solved')
fu.countour_plot_via_mesh(geometry, u, levels = p.contour_levels, PATH = PATH, plot_title = '')

# fu.fenics_plot(u, PATH, plot_title='')

fu.What_time_is_it(t0, "\u03C8(r, z) is plotted")
logger.log_n_output_colored_message(colored_message="'Done'\n", color='red', white_message='')
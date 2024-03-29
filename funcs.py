import MEPHIST_data as MEPH
from imports import *
import time
import math
import point_source_data as psd
import logger
import mshr
import matplotlib.tri as tri
import pylab
from MEPhIST_psi_axis_problem_params import Problem
import MEPhIST_2_problems_problem_params as MEPhIST_exp_profile
# %% Problem parameters
DEFAULT_MESH = 100  # for mesher characterized by 2 params
MESH_DENSITY = 20  # for mesher characterized by 1 param

EPS = 0.05  # when zero maybe inf (1/r)
Z0 = 0.4
R1, Z1 = 0, -Z0  # see Krat's unpublishet article
R2, Z2 = 0.55, Z0
# %% Some consts
DIRICHLET_BOUNDARY = 'DIRICHLET_BOUNDARY'
NEUMANN_BOUNDARY = 'NEUMANN_BOUNDARY'

SOLVE_PLASMA_POINT_SOURCES = 1
SOLVE_PLASMA = 2
SOLVE_POINT_SOURCES = 3
SOLVE_PLASMA_POINT_SOURCES_LINEAR_SOLVER = 4

SOLVE_DICT = {
    1: "SOLVE_PLASMA_POINT_SOURCES",
    2: "SOLVE_PLASMA",
    3: "SOLVE_POINT_SOURCES",
    4: "SOLVE_PLASMA_POINT_SOURCES_LINEAR_SOLVER"
}

M0 = 1.25e-6
# %% Plot stuff
DPI = 200  # quality of plots
TEXT_FILE_U_MAX = "Text_data/func_max"
TEXT_FILE_ERROR = "Text_data/error"
TEXT_FILE_2D_PLOT = "Text_data/2Dplot"
TWOD_PLOT_SAVE_PATH = 'Figures/Post_analyt'
# %% Square influence  stuff
SQ_MIN = 1
SQ_MAX = 9
SQUARE_SIZE_ARRAY = numpy.linspace(
    SQ_MIN, SQ_MAX, 1+int((SQ_MAX-SQ_MIN)/abs(SQ_MIN)))
# %% Fonts for plots
FONT = {'family': "Times New Roman",
        'size': 15}
matplt.rc('font', **FONT)
# %% Funcs


def Form_f_text(A1, A2):
    # A1 = mo*p', A2 = FF'
    # deriviation are calculated using 'sympy' library
    x = sympy.symbols('x[0]')  # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("Right-hand equation side (f): \n", 'magenta') + f_text)

    return f_text


def Twod_plot(psi, x0, y1, y2, path, square_size):
    # y1, y2 - min and max points in an interval of interest,
    # x0 - point along which 2d graph is plotted
    # psi.set_allow_extrapolation(True)
    # avoid hitting points outside the domain
    tol, point_num = 0.001, DEFAULT_MESH + 1
    y = numpy.linspace(y1 + tol, y2 - tol, point_num)

    points = [(x0, y_) for y_ in y]  # create 2D points
    psi_line = numpy.array([psi(point) for point in points])

    Save_2D_data(square_size, numpy.array([y, psi_line]).transpose())

    matplt.plot(y, psi_line, 'k', linewidth=2)  # magnify w
    matplt.grid(True)
    matplt.xlabel('r')
    matplt.ylabel('\u03C8')
    matplt.legend(["Сечение r0: %s, размер квадрата: %s" %
                  (x0, square_size), 'Load'], loc='best')

    time_title = Time_name()

    file_path = 'Figures/%s/%s_%s.png' % (path, time_title, x0)
    matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
    matplt.close()  # close created plot

    print(colored("2D plot saved to PATH: %s" % file_path, 'green'))
    return numpy.amax(psi_line)


def Cross_section_x0(x0, y1, y2, psi, tol, point_num):
    # its here for future
    y = numpy.linspace(y1 + tol, y2 - tol, point_num)

    points = [(x0, y_) for y_ in y]  # create 2D points
    psi_line = numpy.array([psi(point) for point in points])

    return points, psi_line


def Save_2D_data(square_size, data):
    # variable names are self explanatory
    file_path = "%s_%s_%s.txt" % (TEXT_FILE_2D_PLOT, DEFAULT_MESH, square_size)
    with open(file_path, 'w') as file:
        for line in data:
            file.write("%s,%s\n" % (line[0], line[1]))

    print(colored("2D cross-section data saved to PATH: %s" % file_path, 'green'))


def Plot_2D_data_together():  # father
    iteration = 0

    for i in SQUARE_SIZE_ARRAY[::2]:
        # for i in SQUARE_SIZE_ARRAY:
        # variable names are self explanatory
        file_path = "%s_%s_%s.txt" % (TEXT_FILE_2D_PLOT, DEFAULT_MESH, i)
        with open(file_path, "r") as file:  # change to Read_from_file func
            data = [[float(num) for num in line.split(',')] for line in file]

        x = numpy.array(get_column(data, 0))
        u_section = numpy.array(get_column(data, 1))

        if iteration == 0:
            u0_section = u_section
            delta_arr = []
            i_array = []
        else:
            u_section = Level_arrays(u0_section, u_section)

            delta_arr.append((100*(u_section-u0_section)/u0_section).tolist())
            i_array.append(i)

        matplt.plot(x, u_section, linewidth=1,
                    label="Размер квадрата = %s" % int(i))

        iteration = iteration + 1

    matplt.xlabel('r')
    matplt.ylabel('\u03C8')
    matplt.grid(True)
    matplt.legend()

    file_path = "%s/2D_plots_together_%s.png" % (
        TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi=2*DPI, bbox_inches="tight")
    matplt.close()  # close created plot
    print(colored("2D plots saved together at PATH: %s" % file_path, 'green'))

    # now build errors # maybe change it for the better next time!
    iteration = 0
    n = 20
    for i in i_array:
        matplt.plot(x[n:len(x)-n], delta_arr[iteration][n:len(x)-n],
                    linewidth=0.5, label="Размер квадрата = %s" % int(i))
        iteration = iteration + 1
    matplt.xlabel('r')
    matplt.ylabel("\u0394, %")
    matplt.grid(True)
    matplt.legend()

    file_path = "%s/2D_plots_together_error_%s.png" % (
        TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi=2*DPI, bbox_inches="tight")
    matplt.close()  # close created plot

    print(colored("2D plots of error saved together at PATH: %s" % file_path, 'green'))


def Level_arrays(u0, u1):
    # use only in point source boundary influence studies only!
    diff = max(u0) - max(u1)
    u1 = numpy.array(u1) + diff*numpy.ones(len(u1))

    return u1


def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  # get current time to make figure name unique
    return time_title


def Save_figure(f_expr, mesh_r, mesh_z, addition, PATH, plot_title):
    # move to funcs, add missing args, fix save path
    # Plot solution and mesh. Save plot
    # nothing passed to function, because variables are global
    mesh_title = "%sx%s mesh" % (str(mesh_r), str(mesh_z))

    time_title = Time_name()

    # file path+unique time name
    path_my_file = 'Figures/%s/%s' % (PATH, time_title)

    matplt.title(plot_title)  # titled figure for my self

    file_path = "%s.png" % path_my_file
    logger.info(file_path)
    # no title figure for reports
    matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
    matplt.close()  # close created plot
    time.sleep(2)

    print_colored_n_white(
        colored_text="3D countour plot saved to PATH: ", color='green', white_text=file_path)


# def write_data_to_file(folder_name, file_name, x, y, z=[]):
#     file_path = "%s/%s.txt" % (folder_name, file_name)
#     file = open(file_path, "r")  # append write to file mode
    
#     text = "%s,%s,%s\n" % (x, y, z)
#     file.write(text)
#     file.close()

#     print(colored("Data saved to PATH: %s" % file_path, 'green'))


def Write2file_umax_vs_square_size(mesh_r, mesh_z, u_max, default_mesh_size):
    file_path = "%s_vs_square_mesh_%s.txt" % (
        TEXT_FILE_U_MAX, default_mesh_size)
    file = open(file_path, "a")  # append write to file mode

    text = "%s,%s,%s\n" % (mesh_r, mesh_z, u_max)
    file.write(text)
    file.close()

    print(colored("Data saved to PATH: %s" % file_path, 'green'))


def get_column(matrix, col):
    return [row[col] for row in matrix]


# u max as a function of mesh parameters on the same solution area
def Plot_umax_vs_def_mesh(name):
    with open("%s.txt" % TEXT_FILE_U_MAX, "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]

    mesh = get_column(data, 0)
    u_max = get_column(data, 2)

    matplt.scatter(mesh, u_max, linewidth=2)  # magnify w
    # matplt.legend(["u_max vs default mesh size"], loc='best')
    matplt.grid(True)
    matplt.xlabel('Размер сетки')
    matplt.ylabel('\u03C8')

    matplt.savefig("Figures/umax_vs_mesh_%s.png" %
                   name, dpi=DPI, bbox_inches="tight")

    matplt.close()  # close created plot


# u max as a function of solution square size
def Plot_umax_vs_square_size(name, default_mesh_size):
    with open("%s_vs_square_mesh_%s.txt" % (TEXT_FILE_U_MAX, default_mesh_size), "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]

    mesh = get_column(data, 0)
    u_max = get_column(data, 2)

    matplt.scatter(mesh, u_max, linewidth=2)
    # matplt.legend(["u_max vs solution square size"], loc='best')
    matplt.grid(True)
    matplt.xlabel("Размер области \u03A9")
    matplt.ylabel('\u03C8')
    # matplt.title("Default mesh size: %d" % (default_mesh_size)) # titled figure for my self

    matplt.savefig("Figures/umax_vs_sq_sz_%s_%s.png" %
                   (name, default_mesh_size), dpi=DPI, bbox_inches="tight")

    matplt.close()  # close created plot


def What_time_is_it(t0, message):
    print(colored("\tTime elapsed = %f (%s)" %
          (time.time() - t0, message), 'blue'))


def Analyt_sol(c, A1, A2):
    x = sympy.symbols('x[0]')  # r coordinate
    z = sympy.symbols('x[1]')  # r coordinate
    # sympy.log
    psi_p = A1 * pow(x, 4) + A2 * pow(z, 2)  # private solution
    psi_gen = \
        c[0] + \
        c[1] * pow(x, 2) + \
        c[2] * (pow(x, 4) - 4*pow(x, 2)*pow(z, 2)) + \
        c[3] * (-pow(z, 2))  # general solution the rest of the 4th term is defined in MyLog(c) func
    # c[3] * (pow(x, 2)*sympy.log(x)- pow(z, 2)) # general solution
    #pow(x, 2)*sympy.log(x)

    my_log = MyLog(c)

    psi_text = sympy.printing.ccode(psi_p + psi_gen)
    psi_p_text = sympy.printing.ccode(psi_p)

    # final_sol = psi_text
    final_sol = psi_text + ' + ' + my_log
    print(colored("Private solution: \n", 'magenta') + psi_p_text)
    print(colored("Analytical solution: \n", 'magenta') + final_sol)
    # print(colored("Analytical solution: \n", 'magenta') + psi_text)

    return final_sol


def MyLog(c):
    x = sympy.symbols('x[0]')  # r coordinate
    pre_log = c[3] * pow(x, 2)

    pre_log_text = sympy.printing.ccode(pre_log)
    # assemble function of the point source
    log_text = "%s*std::log(%s)" % (pre_log_text, 'x[0]')
    print(colored("Problem term in analyt solution: \n", 'magenta') + log_text)
    # c[3] * (pow(x, 2)*log(x)) # general solution

    return log_text


def CreatePointSource(r, I, disp):
    x = sympy.symbols('x[0]')  # r coordinate
    z = sympy.symbols('x[1]')  # r coordinate

    # in sympy write stuff that works
    pre_exp = 2 * M0/pi/math.pow(disp, 2)/math.erfc(-r[0]/disp) * I * x
    inner_exp = - (pow(x - r[0], 2) + pow(z - r[1], 2)) / \
        math.pow(disp, 2)  # in sympy write stuff that works
    pre_exp_text = sympy.printing.ccode(pre_exp)  # transfer it to text
    inner_exp_text = sympy.printing.ccode(inner_exp)  # transfer it to text

    # assemble function of the point source
    point_source_text = "%s*std::exp(%s)" % (pre_exp_text, inner_exp_text)
    logger.log_n_output(point_source_text, 'white')
    # point_source_text = point_source_text.replace('pow', 'std::pow') # reason being faulty fenics namespace

    return point_source_text


def ArrayOfPointSources(pnt_src_data):
    # create an array of all point source text expressions

    logger.info('Point sources params:')
    logger.info('r')
    logger.info(pnt_src_data.r)
    logger.info('I & disp')
    logger.info(pnt_src_data.i_disp)
    logger.log_n_output("Point sources:", 'magenta')

    pnt_src_text = []
    for i in range(len(pnt_src_data.r)):
        pnt_src_text.append(CreatePointSource(
            pnt_src_data.r[i], pnt_src_data.i_disp[i][0], pnt_src_data.alpha * pnt_src_data.i_disp[i][1]))

    return pnt_src_text


def Array_Expression(text_array):
    expression_array = [None]*len(text_array)
    for i in range(len(text_array)):
        expression_array[i] = Expression(text_array[i], degree=2)

    return expression_array

def Mesh_to_xml():
    file = File('Mesh/file.xml')

# def D_config(smoothness):
#     tol = 0.0001
#     t = numpy.linspace(tol, math.pi - tol, smoothness)
#     x = numpy.ones(smoothness) + EPS*numpy.cos(t + ALPHA*numpy.sin(t))
#     z = EPS*KAPPA*numpy.sin(t)

#     x = numpy.append(x, numpy.flip(x)) # move it along x axis!
#     z = numpy.append(z, numpy.flip(-z))
#     matplt.plot(x, z)
#     matplt.grid()
#     matplt.show()


def Write2file_errors(mesh_r, mesh_z, err_L2, err_max):
    file_path = "%svs%s.txt" % (TEXT_FILE_ERROR, 'mesh')
    file = open(file_path, "a")  # append write to file mode

    text = "%s,%s,%s\n" % (mesh_r, err_max, err_L2)
    file.write(text)
    file.close()

    print(colored("Data saved to PATH: %s" % file_path, 'green'))


# u max as a function of mesh parameters on the same solution area
def Plot_error_vs_mesh(name):
    with open("%svsmesh.txt" % TEXT_FILE_ERROR, "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]

    mesh = get_column(data, 0)
    err_max = get_column(data, 1)
    err_L2 = get_column(data, 2)

    matplt.semilogy(mesh, err_max, 'o', linewidth=2)  # magnify w
    # matplt.legend(["u_max vs default mesh size"], loc='best')
    matplt.grid(True)
    matplt.xlabel('Размер сетки')
    matplt.ylabel('Максимальная ошибка')

    file_path = 'Figures/Post_analyt/analt_errorvsmesh.png'
    matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
    print(colored("2D plot saved to PATH: %s" % file_path, 'green'))

    matplt.close()  # close created plot


def To_float(arr):
    arr_str = []
    for i in arr:
        arr_str.append(str(i))

    return arr_str


def Hand_input(p_pow, F_pow):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')

    # flux function #think tomorrow how to define argument psi!
    psi = sympy.symbols('u')
    # r coordinate. used for easy writing of expressions
    x = sympy.symbols('x[0]')

    p_psi = pow(psi/M.psi_axis, int(p_pow))  # pressure function
    F_psi_2 = 1 - pow(psi/M.psi_axis, int(F_pow))  # poloidal current function

    dp_psi = sympy.diff(p_psi, psi)  # pressure and F deriviation
    dF_psi_2 = sympy.diff(F_psi_2, psi)  # compiler breaks when

    f_text = (M0 * pow(x, 2) * M.p_axis * dp_psi + 0.5 *
              M.F0_2 * dF_psi_2)  # right hand expression

    f_text = sympy.printing.ccode(f_text)
    # reason being faulty fenics namespace
    f_text = f_text.replace('exp', 'std::exp')
    # f_text = f_text.replace('pow', 'std::pow') # reason being faulty fenics namespace

    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output(f_text, 'white')

    return f_text


def Initial_guess_for_u(u, const):
    for i in range(len(u.vector())):
        u.vector()[i] = float(const)
    logger.log_n_output_colored_message(
        colored_message="Initial guess for u: ", color='green', white_message=str(const))
    return u


def level_u(u, const):
    for i in range(len(u.vector())):
        u.vector()[i] = float(const) + u.vector()[i]
    logger.log_n_output_colored_message(
        colored_message="Boundary condition added for u: ", color='green', white_message=str(const))
    return u


def Neumann_boundary(x, on_boundary):
    tol = 1e-10
    return on_boundary and (abs(x[1] - Z1) < tol or abs(x[1] - Z2) < tol or abs(x[0] - R2) < tol)


def Dirichlet_boundary(x, on_boundary):
    return on_boundary


def print_colored(text, color):
    print(colored(text, color))


def refine_mesh(mesh, domains_amount):
    cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
    cell_markers.set_all(False)
    for i in range(domains_amount):
        submesh = SubMesh(mesh, i+1)
        [cell_markers, global_mesh_index] = refine_subdomain(
            mesh, submesh, cell_markers)  # domain numbering starts wwith 1

    mesh = refine(mesh, cell_markers)
    logger.info("Mesh refined: %s" % str(global_mesh_index))
    logger.info("After refinement. Number of cells: %d, Number of vertices: %d" % (
        mesh.num_cells(), mesh.num_vertices()))

    return mesh


def refine_subdomain(mesh, submesh, cell_markers):
    bound_box = mesh.bounding_box_tree()
    global_mesh_index = []

    for cell in cells(submesh):
        global_mesh_index.append(
            bound_box.compute_first_entity_collision(cell.midpoint()))

    for cell in cells(mesh):
        if cell.index() in global_mesh_index:
            cell_markers[cell] = True
        else:
            cell_markers[cell] = False

    return cell_markers, global_mesh_index


def Create_Subdomain(r, disp, segments):
    circle = mshr.Circle(Point(r[0], r[1]), disp*3, segments=segments)
    logger.info("Created subdomain: %s, r=%e, segments=%d" %
                (str(r), disp*3, segments))

    return circle


def Create_Subdomains(alpha, segments):
    circle_list = []
    ps_d = psd.PointSource(alpha)

    for i in range(len(ps_d.r)):
        if abs(ps_d.i_disp[i][0]) > 1e-6:
            circle_list.append(Create_Subdomain(
                ps_d.r[i], ps_d.i_disp[i][1] * ps_d.alpha, segments))

    return circle_list


def Set_Subdomains(domain, alpha, segments):
    circle_list = Create_Subdomains(alpha, segments)
    for i in range(len(circle_list)):
        domain.set_subdomain(i+1, circle_list[i])

    return domain, len(circle_list)


def plot_mesh(mesh, path):
    plot(mesh)

    matplt.xlabel("r, м")
    matplt.ylabel("z, м")
    matplt.xticks([0.1, 0.2, 0.3, 0.4, 0.5])
    matplt.gca().set_aspect("equal")

    Save_figure('', 100, 100, '', path, "")

    return 0


def plasma_sources_coefficients_pow_2(p_correction, F_correction):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')

    p_coeff = 2 * M0 * M.p_axis / M.psi_axis**2 * p_correction
    F_2_coeff = -M.F0_2 / M.psi_axis**2 * F_correction

    logger.log_n_output_colored_message(
        colored_message="p_correction = ", color='green', white_message=str(p_correction))
    logger.log_n_output_colored_message(
        colored_message="F_2_correction = ", color='green', white_message=str(F_correction))

    logger.log_n_output_colored_message(
        colored_message="p_coeff = ", color='green', white_message=str(p_coeff))
    logger.log_n_output_colored_message(
        colored_message="F_2_coeff = ", color='green', white_message=str(F_2_coeff))

    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output("%s*pow(x[0], 2)*u + %s*u" %
                        (p_coeff, F_2_coeff), 'white')

    return p_coeff, F_2_coeff


def plasma_sources_coefficients_pow_2_iteration(p_correction, F_correction, psi_axis):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')

    p_coeff = 2 * M0 * M.p_axis / psi_axis**2 * p_correction
    F_2_coeff = -M.F0_2 / psi_axis**2 * F_correction

    logger.log_n_output_colored_message(
        colored_message="p_correction = ", color='green', white_message=str(p_correction))
    logger.log_n_output_colored_message(
        colored_message="F_2_correction = ", color='green', white_message=str(F_correction))

    logger.log_n_output_colored_message(
        colored_message="p_coeff = ", color='green', white_message=str(p_coeff))
    logger.log_n_output_colored_message(
        colored_message="F_2_coeff = ", color='green', white_message=str(F_2_coeff))

    logger.log_n_output_colored_message(
        colored_message="psi axis (by me) = ", color='green', white_message=str(psi_axis))

    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output("%s*pow(x[0], 2)*u + %s*u" %
                        (p_coeff, F_2_coeff), 'white')

    return p_coeff, F_2_coeff


def plasma_sources_coefficients_exp_profile(p_correction, F_correction, psi_correction):
    M = MEPH.MEPhIST()
    exp_prof = MEPhIST_exp_profile.Problem()

    psi1 = exp_prof.psi_axis
    psi2 = exp_prof.psi_pl_edge

    p_0 = M.p_axis
    F0_2 = M.F0_2

    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')

    p_coeff = M0 * p_0 * p_correction / (psi1 - psi2) / psi_correction
    F_2_coeff = -F0_2 * F_correction / (psi1 - psi2) / psi_correction

    logger.log_n_output_colored_message(
        colored_message="p_correction = ", color='green', white_message=str(p_correction))
    logger.log_n_output_colored_message(
        colored_message="F_2_correction = ", color='green', white_message=str(F_correction))

    logger.log_n_output_colored_message(
        colored_message="psi1 = ", color='green', white_message=str(psi1))
    logger.log_n_output_colored_message(
        colored_message="psi2 = ", color='green', white_message=str(psi2))

    logger.log_n_output_colored_message(
        colored_message="p_coeff = ", color='green', white_message=str(p_coeff))
    logger.log_n_output_colored_message(
        colored_message="F_2_coeff = ", color='green', white_message=str(F_2_coeff))

    logger.log_n_output_colored_message(
        colored_message="psi axis (by me) = ", color='green', white_message=str(psi_correction))

    return p_coeff, F_2_coeff


def spheromak_point(r, R, alpha):
    z = 1/2/alpha * math.sqrt(2 * R**2 - (r - R)**2)
    return z


def plot_spheromak_boundary(R, alpha, smoothness):
    r_array = numpy.linspace((1 + math.sqrt(2)) * R,
                             (1 - math.sqrt(2)) * R, smoothness)
    z_array = []
    for r in r_array:
        z_array.append(spheromak_point(r, R, alpha))

    z_array = numpy.array(z_array)

    r_array = numpy.append(r_array, numpy.flip(r_array))  # замкнуть кривую
    z_array = numpy.append(z_array, numpy.flip(-z_array))  # замкнуть кривую
    return r_array, z_array


def spheromak_boundary(R, alpha, smoothness):
    [r_array, z_array] = plot_spheromak_boundary(R, alpha, smoothness)
    boundary_geometry = []
    for i in range(len(r_array)):
        boundary_geometry.append(Point(r_array[i], z_array[i]))

    return mshr.Polygon(boundary_geometry)


def spheromak_pressure(psi_0, R, alpha):
    L = Expression(
        "pow(x[0], 2) / pow(%s, 4) * (1 + pow(%s, 2)) * %s" % (R, alpha, psi_0), degree=2)
    logger.log_n_output("Right hand part:", 'red')
    logger.log_n_output(L._cppcode, 'white')
    return L


def countour_plot_via_mesh(geometry, u, levels, PATH,
                           plot_title='',
                           current_disp=[0, 0],
                           do_plasma_centre = False,
                           plt_vessel=False,
                           xticks_array = [],
                           yticks_array = [],
                           grid = False,
                           colorbar = False):

    u_min = u.vector()[:].min()
    u_max = u.vector()[:].max()
    
    if u_min == u_max:
        logger.log_n_output(message="Trivial solution. u = %s" %
                            u_max, color='red')
    else:
        triang = tri.Triangulation(
            *geometry.mesh.coordinates().reshape((-1, 2)).T, triangles=geometry.mesh.cells())
        u_array = u.compute_vertex_values(geometry.mesh)

        if do_plasma_centre == True:        
            matplt.scatter(current_disp[0], current_disp[1], c='r', linewidth=2.5, zorder=3) if current_disp != 0 else 0
            
        fig = matplt.tricontour(triang, u_array, levels, zorder=2)
        matplt.gca().set_aspect("equal")
        
        matplt.xlim(geometry.plot_domain[0], geometry.plot_domain[1])
        matplt.ylim(geometry.plot_domain[2], geometry.plot_domain[3])
        matplt.xlabel("r, м")
        matplt.ylabel("z, м")
        
        if plt_vessel == True:
            matplt.plot(geometry.outer_vessel_contour[0],
                        geometry.outer_vessel_contour[1],
                        c='k', linewidth=1, zorder=1)
            matplt.plot(geometry.inner_vessel_contour[0],
                        geometry.inner_vessel_contour[1],
                        c='k', linewidth=1, zorder=1)
            
        if xticks_array != []:
            matplt.xticks(numpy.array(xticks_array))
        if yticks_array != []:
            matplt.xticks(numpy.array(yticks_array))

        if grid == True:
            matplt.grid(True)
        
        if colorbar == True:
            matplt.colorbar(fig).set_label("\u03C8(r, z), Вб")

        logger.log_n_output_colored_message(
            colored_message="u_max = ", color='green', white_message=str(u_max))
        logger.log_n_output_colored_message(
            colored_message="u_min = ", color='green', white_message=str(u_min))
        logger.log_n_output_colored_message(
            colored_message="u_max-u_min = ", color='green', white_message=str(u_max-u_min))
        logger.info("levels = %s" % str(levels))

        save_contour_plot(PATH, plot_title)

        return u_max
    
# def countour_plot_via_mesh_nocolorbar(geometry, u, levels, PATH,
#                                       plot_title='', current_disp=0,
#                                       plt_vessel=False):
#     u_min = u.vector()[:].min()
#     u_max = u.vector()[:].max()
#     if u_min == u_max:
#         logger.log_n_output(message="Trivial solution. u = %s" %
#                             u_max, color='red')
#     else:
#         triang = tri.Triangulation(
#             *geometry.mesh.coordinates().reshape((-1, 2)).T, triangles=geometry.mesh.cells())
#         u_array = u.compute_vertex_values(geometry.mesh)

#         if plt_vessel == True:
#             matplt.plot(geometry.outer_vessel_contour[0],
#                         geometry.outer_vessel_contour[1],
#                         c='k', linewidth=1)
#             matplt.plot(geometry.inner_vessel_contour[0],
#                         geometry.inner_vessel_contour[1],
#                         c='k', linewidth=1)
        
#         matplt.xticks(numpy.array([0.1, 0.2, 0.3, 0.4, 0.5]))
#         # matplt.xlim([geometry.domain_boundary_coordinates[0]])
#         matplt.grid(True)
#         matplt.tricontour(triang, u_array, levels)
        
#         matplt.scatter(current_disp, 0, c='r', linewidth=2.5) if current_disp != 0 else 0

#         matplt.xlabel("r, м")
#         matplt.ylabel("z, м")
#         matplt.gca().set_aspect("equal")

#         save_contour_plot(PATH, plot_title)

#         return u_max


def fenics_plot(problem, u, PATH,
                plot_title='',
                xticks=[],
                yticks=[],
                limits=0,
                colorbar=False,
                show=False):
    # fig = plot(u, linewidth=0.5)
    fig = plot(u)
    if colorbar == True:
        pylab.colorbar(fig).set_label("\u03C8(r, z), Вб")

    if limits == 0:
        matplt.xlim(problem.domain_geometry0[0], problem.domain_geometry0[1])
        matplt.ylim(problem.domain_geometry0[2], problem.domain_geometry0[3])
    else:
        matplt.xlim(problem.domain_geometry[0], problem.domain_geometry[1])
        matplt.ylim(problem.domain_geometry[2], problem.domain_geometry[3])
    
    if xticks != []:
        matplt.xticks(numpy.array(xticks))
    
    if yticks != []:
        matplt.xticks(numpy.array(yticks))
    
    matplt.gca().set_aspect("equal")
    matplt.xlabel("r, м")
    matplt.ylabel("z, м")
    matplt.grid(True)
    
    if show != False:
        matplt.show()
    save_contour_plot(PATH, plot_title)


def save_contour_plot(PATH, plot_title):
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title)
    file_path = "%s.png" % path_my_file
    
    logger.info(file_path)

    matplt.title(plot_title)
    matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
    matplt.close()

    print_colored_n_white(
        colored_text="3D countour plot saved to PATH: ",
        color='green',
        white_text=file_path)
    time.sleep(1)


def print_colored_n_white(colored_text, color, white_text):
    print(colored(colored_text, color) + white_text)


def ErrorEstimate(u, u_D, mesh):
    error_L2 = errornorm(u_D, u, 'L2')

    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = numpy.max(numpy.abs(vertex_values_u_D - vertex_values_u))

    logger.log_n_output_colored_message(
        colored_message="error_L2 = ", color='red', white_message=str(error_L2))
    logger.log_n_output_colored_message(
        colored_message="error_max = ", color='red', white_message=str(error_max))

    return error_L2, error_max
# %%


def plot_1D(PATH, u, geometry):

    plot(u)

    matplt.grid("True")
    matplt.xlim(geometry.a, geometry.b)
    matplt.xlabel("r, м")
    matplt.ylabel("p(r), Па")

    logger.print_colored_n_white(
        colored_text="u_max = ", color='green', white_text=str(u.vector()[:].max()))
    logger.print_colored_n_white(
        colored_text="u_min = ", color='green', white_text=str(u.vector()[:].min()))
    logger.print_colored_n_white(colored_text="u_max-u_min = ", color='green',
                                 white_text=str(u.vector()[:].max() - u.vector()[:].min()))

    save_contour_plot(PATH, "")


def multiply_u_by_const(u, const):
    for i in range(len(u.vector())):
        u.vector()[i] = float(const)*u.vector()[i]
    logger.log_n_output_colored_message(
        colored_message="u is multiplied by: ", color='green', white_message=str(const))
    return u


def plot_Dina_results(PATH):
    problem = Problem()

    x, y, z = numpy.genfromtxt(r'psi.dat', unpack=True)

    r_lim = [x.min(), x.max()]
    z_lim = [y.min(), y.max()]

    to_m = 1e-2  # to meters units conversion

    levels_min = 0
    levels_max = 0.055  # default good 0.055

    # levels_amount = 12 # default good 3*(1 + int(100*(abs(levels_min)+abs(levels_max))))
    # print_colored_n_white(colored_text="levels values = ", color='green', white_text=str(levels_amount))
    levels = 20  # list(numpy.linspace(levels_min, levels_max, levels_amount))

    print_colored_n_white(colored_text="r min = ",
                          color='green', white_text=str(r_lim[0]))
    print_colored_n_white(colored_text="r max = ",
                          color='green', white_text=str(r_lim[1]))

    print_colored_n_white(colored_text="z min = ",
                          color='green', white_text=str(z_lim[0]))
    print_colored_n_white(colored_text="z max = ",
                          color='green', white_text=str(z_lim[1]))

    print_colored_n_white(colored_text="u min = ",
                          color='green', white_text=str(z.min()))
    print_colored_n_white(colored_text="u max = ",
                          color='green', white_text=str(z.max()))

    matplt.tricontour(to_m*x, to_m*y, z, levels=levels)
    matplt.colorbar().set_label("\u03C8(r, z), Вб")
    # format="%.2f"

    matplt.xlim(problem.domain_geometry[0], problem.domain_geometry[1])
    matplt.ylim(problem.domain_geometry[2], problem.domain_geometry[3])
    matplt.xlabel("r, м")
    matplt.ylabel("z, м")

    matplt.grid(True, which='both')
    matplt.gca().set_aspect("equal")

    save_contour_plot(PATH, "")


def plot_error_vs_mesh_density(mesh_array, errors, PATH):
    matplt.semilogy(mesh_array, errors, 'o', linewidth=2)  # magnify w
    matplt.grid(True)
    matplt.xlabel('Плотность сетки')
    matplt.ylabel('Максимальная ошибка, Вб')

    save_contour_plot(PATH, "")


def save_errors_to_file(mesh_array, errors, file_name):
    time_title = Time_name()
    file_path = "Errors/%s_%s.txt" % (file_name, time_title)
    with open(file_path, 'w') as file:
        for i in range(len(mesh_array)):
            file.write("%s,%s\n" % (mesh_array[i], errors[i]))

    print_colored_n_white(colored_text="Data saved to PATH: ",
                          color='green', white_text=file_path)


# u max as a function of mesh parameters on the same solution area
def plot_error_vs_mesh_from_file(folder_name, file_name, x_lim, PATH):
    with open("%s/%s.txt" % (folder_name, file_name), "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]

    mesh = get_column(data, 0)
    error = get_column(data, 1)
    
    # curve_fit(approximate_errors, mesh, error)    

    matplt.loglog(mesh, error, 'o', linewidth=2)  # magnify w
    matplt.xlim(x_lim[0], x_lim[1])
    matplt.grid(True)
    matplt.xlabel('Плотность сетки')
    matplt.ylabel('Максимальная ошибка, Вб')

    save_contour_plot(PATH, '')
    
def approximate_errors(x, a, b, c):
    return a * numpy.exp(-b * x) + c

def shrink_contour(folder_name, file_name, shrunk_to_point, alpha_x, alpha_z):
    [x, z] = read_from_file(folder_name, file_name)
    
    x_shrinking = shrunk_to_point[0] * numpy.ones(len(x))
    z_shrinking = shrunk_to_point[1] * numpy.ones(len(x))
    
    # print("%s" % [x_mid[0], z_mid[0]])
    x_shrunk = alpha_x*x + x_shrinking*(1-alpha_x)
    z_shrunk = alpha_z*z + z_shrinking*(1-alpha_z)
    
    matplt.plot(x, z)
    matplt.plot(x_shrunk, z_shrunk)
    matplt.xticks(numpy.linspace(0.1, 0.4, 7))
    matplt.grid(True)
    matplt.gca().set_aspect("equal")
    # matplt.show()
    
    from geometry import Geometry
    
    ge = Geometry()
    data = numpy.transpose([x_shrunk,z_shrunk])
    ge.write_data_to_file('Data', 'MEPHIST_vessel_inner_surface', data)    
        
def read_from_file(folder_name, file_name):
    file_path = "%s/%s.txt" % (folder_name, file_name)
    with open(file_path, "r") as file: # change to Read_from_file func
        data = [[float(num) for num in line.split(',')] for line in file]
    
    x = numpy.array(get_column(data, 0))
    z = numpy.array(get_column(data, 1))

    return x, z    

def write_data_to_file(folder_name, file_name, data):
    if len(data) == 2:
        data = numpy.transpose(data)
    file_path = "%s/%s.txt" % (folder_name, file_name)
    with open(file_path, 'w') as file:
        for line in data:
            file.write("%s,%s\n" % (line[0], line[1]))
from matplotlib.pyplot import text
import MEPHIST_data as MEPH
from imports import *
import time
import math
import point_source_data as psd
import logger
import mshr
import matplotlib.tri as tri
#%% Problem parameters
DEFAULT_MESH = 100 # for mesher characterized by 2 params
MESH_DENSITY = 20 # for mesher characterized by 1 param

EPS = 0.05 # when zero maybe inf (1/r)
Z0 = 0.4
R1, Z1 = 0, -Z0 # see Krat's unpublishet article
R2, Z2 = 0.55, Z0
#%% Some consts
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
#%% Plot stuff
DPI = 200 # quality of plots
TEXT_FILE_U_MAX = "Text_data/func_max"
TEXT_FILE_ERROR = "Text_data/error"
TEXT_FILE_2D_PLOT = "Text_data/2Dplot"
TWOD_PLOT_SAVE_PATH = 'Figures/Post_analyt'
#%% Square influence  stuff
SQ_MIN = 1
SQ_MAX = 9
SQUARE_SIZE_ARRAY = numpy.linspace(SQ_MIN, SQ_MAX, 1+int((SQ_MAX-SQ_MIN)/abs(SQ_MIN)))
#%% Fonts for plots
FONT = {'family' : "Times New Roman",
        'size' : 15}
matplt.rc('font', **FONT)
#%% Funcs
def Form_f_text(A1, A2):
    #A1 = mo*p', A2 = FF'
    #deriviation are calculated using 'sympy' library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("Right-hand equation side (f): \n", 'magenta') + f_text)

    return f_text

def Twod_plot(psi, x0, y1, y2, path, square_size): 
    # y1, y2 - min and max points in an interval of interest, 
    # x0 - point along which 2d graph is plotted
    # psi.set_allow_extrapolation(True)
    tol, point_num = 0.001, DEFAULT_MESH + 1  # avoid hitting points outside the domain
    y = numpy.linspace(y1 + tol, y2 - tol, point_num)
    
    points = [(x0, y_) for y_ in y]  # create 2D points
    psi_line = numpy.array([psi(point) for point in points])
    
    Save_2D_data(square_size, numpy.array([y, psi_line]).transpose())
    
    matplt.plot(y, psi_line, 'k', linewidth=2)  # magnify w
    matplt.grid(True)
    matplt.xlabel('r')
    matplt.ylabel('\u03C8')
    matplt.legend(["Сечение r0: %s, размер квадрата: %s" % (x0, square_size), 'Load'], loc='best')
    
    time_title = Time_name()
    
    file_path = 'Figures/%s/%s_%s.png' % (path, time_title, x0)
    matplt.savefig(file_path, dpi = DPI, bbox_inches="tight")
    matplt.close() # close created plot
    
    print(colored("2D plot saved to PATH: %s" % file_path, 'green'))
    return numpy.amax(psi_line)

def Cross_section_x0(x0, y1, y2, psi, tol, point_num):
    # its here for future
    y = numpy.linspace(y1 + tol, y2 - tol, point_num)
    
    points = [(x0, y_) for y_ in y]  # create 2D points
    psi_line = numpy.array([psi(point) for point in points])
    
    return points, psi_line
   
def Save_2D_data(square_size, data):
    file_path = "%s_%s_%s.txt" % (TEXT_FILE_2D_PLOT, DEFAULT_MESH, square_size) # variable names are self explanatory
    with open(file_path,'w') as file:
        for line in data:
            file.write("%s,%s\n" % (line[0], line[1]))
    
    print(colored("2D cross-section data saved to PATH: %s" % file_path, 'green'))
    
def Plot_2D_data_together(): #father
    iteration = 0
    
    for i in SQUARE_SIZE_ARRAY[::2]:
    # for i in SQUARE_SIZE_ARRAY:
        file_path = "%s_%s_%s.txt" % (TEXT_FILE_2D_PLOT, DEFAULT_MESH, i) # variable names are self explanatory
        with open(file_path, "r") as file: # change to Read_from_file func
            data = [[float(num) for num in line.split(',')] for line in file]
        
        x = numpy.array(Column(data, 0))
        u_section = numpy.array(Column(data, 1))
        
        
        if iteration == 0:
            u0_section = u_section
            delta_arr = []
            i_array = []
        else:
            u_section = Level_arrays(u0_section, u_section)
            
            delta_arr.append((100*(u_section-u0_section)/u0_section).tolist())
            i_array.append(i)
        
        
        matplt.plot(x, u_section, linewidth=1, label="Размер квадрата = %s" % int(i))
        
        iteration = iteration + 1
    
    
    matplt.xlabel('r')
    matplt.ylabel('\u03C8')
    matplt.grid(True)
    matplt.legend()

    file_path = "%s/2D_plots_together_%s.png" % (TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi = 2*DPI, bbox_inches="tight")
    matplt.close() # close created plot
    print(colored("2D plots saved together at PATH: %s" % file_path, 'green'))
    
    # now build errors # maybe change it for the better next time!
    iteration = 0
    n = 20
    for i in i_array: 
        matplt.plot(x[n:len(x)-n], delta_arr[iteration][n:len(x)-n], linewidth=0.5, label="Размер квадрата = %s" % int(i))
        iteration = iteration + 1
    matplt.xlabel('r')
    matplt.ylabel("\u0394, %")
    matplt.grid(True)
    matplt.legend()
    
    file_path = "%s/2D_plots_together_error_%s.png" % (TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi = 2*DPI, bbox_inches="tight")
    matplt.close() # close created plot
    
    print(colored("2D plots of error saved together at PATH: %s" % file_path, 'green'))

def Level_arrays(u0, u1):
    # use only in point source boundary influence studies only!
    diff = max(u0) - max(u1)
    u1 = numpy.array(u1) + diff*numpy.ones(len(u1))
    
    return u1
    
def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  #get current time to make figure name unique
    return time_title

def Save_figure(f_expr, mesh_r, mesh_z, addition, PATH, plot_title):
    # move to funcs, add missing args, fix save path 
    # Plot solution and mesh. Save plot
    #nothing passed to function, because variables are global
    mesh_title = "%sx%s mesh" % (str(mesh_r), str(mesh_z))
    
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title) # file path+unique time name

    matplt.title(plot_title) # titled figure for my self
    
    file_path = "%s.png" % path_my_file
    logger.info(file_path)
    matplt.savefig(file_path, dpi = DPI, bbox_inches="tight") #no title figure for reports
    matplt.close() # close created plot
    time.sleep(2)
    
    print_colored_n_white(colored_text="3D countour plot saved to PATH: ", color='green', white_text=file_path)
    
def Write2file_umax_vs_def_mesh(mesh_r, mesh_z, u_max):
    file_path = "%s.txt" % TEXT_FILE_U_MAX
    file = open(file_path, "a") # append write to file mode
    
    text = "%s,%s,%s\n" % (mesh_r, mesh_z, u_max)
    file.write(text)
    file.close()
    
    print(colored("Data saved to PATH: %s" % file_path, 'green'))

def Write2file_umax_vs_square_size(mesh_r, mesh_z, u_max, default_mesh_size):
    file_path = "%s_vs_square_mesh_%s.txt" % (TEXT_FILE_U_MAX, default_mesh_size)
    file = open(file_path, "a") # append write to file mode
    
    text = "%s,%s,%s\n" % (mesh_r, mesh_z, u_max)
    file.write(text)
    file.close()
    
    print(colored("Data saved to PATH: %s" % file_path, 'green'))
    
    
def Column(matrix, col):
    return [row[col] for row in matrix]
    
def Plot_umax_vs_def_mesh(name): # u max as a function of mesh parameters on the same solution area
    with open("%s.txt" % TEXT_FILE_U_MAX, "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]
        
    mesh = Column(data, 0) 
    u_max = Column(data, 2)
    
    matplt.scatter(mesh, u_max, linewidth=2)  # magnify w
    # matplt.legend(["u_max vs default mesh size"], loc='best')
    matplt.grid(True)
    matplt.xlabel('Размер сетки')
    matplt.ylabel('\u03C8')
    
    matplt.savefig("Figures/umax_vs_mesh_%s.png" % name, dpi = DPI, bbox_inches="tight")
    
    matplt.close() # close created plot

def Plot_umax_vs_square_size(name, default_mesh_size): # u max as a function of solution square size
    with open("%s_vs_square_mesh_%s.txt" % (TEXT_FILE_U_MAX, default_mesh_size), "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]
        
    mesh = Column(data, 0) 
    u_max = Column(data, 2)
    
    matplt.scatter(mesh, u_max, linewidth=2)
    # matplt.legend(["u_max vs solution square size"], loc='best')
    matplt.grid(True)
    matplt.xlabel("Размер области \u03A9")
    matplt.ylabel('\u03C8')
    # matplt.title("Default mesh size: %d" % (default_mesh_size)) # titled figure for my self
    
    matplt.savefig("Figures/umax_vs_sq_sz_%s_%s.png" % (name, default_mesh_size), dpi = DPI, bbox_inches="tight")
    
    matplt.close() # close created plot
    
def What_time_is_it(t0, message):
    print(colored("\tTime elapsed = %f (%s)" % (time.time() - t0, message), 'blue'))
    
def Analyt_sol(c, A1, A2):
    x = sympy.symbols('x[0]') # r coordinate
    z = sympy.symbols('x[1]') # r coordinate
    #sympy.log
    psi_p = A1 * pow(x, 4) + A2 * pow(z, 2) #private solution
    psi_gen = \
        c[0] + \
        c[1] * pow(x, 2) + \
        c[2] * (pow(x, 4) - 4*pow(x, 2)*pow(z, 2)) + \
        c[3] * (-pow(z, 2)) # general solution the rest of the 4th term is defined in MyLog(c) func
        #c[3] * (pow(x, 2)*sympy.log(x)- pow(z, 2)) # general solution
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
    x = sympy.symbols('x[0]') # r coordinate
    pre_log = c[3] * pow(x, 2)
    
    pre_log_text = sympy.printing.ccode(pre_log)
    log_text = "%s*std::log(%s)" % (pre_log_text, 'x[0]') # assemble function of the point source
    print(colored("Problem term in analyt solution: \n", 'magenta') + log_text)
    #c[3] * (pow(x, 2)*log(x)) # general solution
    
    return log_text

def CreatePointSource(r, I, disp):
    x = sympy.symbols('x[0]') # r coordinate
    z = sympy.symbols('x[1]') # r coordinate

    pre_exp = 2 * M0/pi/math.pow(disp, 2)/math.erfc(-r[0]/disp) * I * x # in sympy write stuff that works
    inner_exp = - (pow(x - r[0], 2) + pow(z - r[1], 2)) / math.pow(disp, 2) # in sympy write stuff that works
    pre_exp_text = sympy.printing.ccode(pre_exp) # transfer it to text
    inner_exp_text = sympy.printing.ccode(inner_exp) # transfer it to text
    
    point_source_text = "%s*std::exp(%s)" % (pre_exp_text, inner_exp_text) # assemble function of the point source
    logger.log_n_output(point_source_text,'white')
    # point_source_text = point_source_text.replace('pow', 'std::pow') # reason being faulty fenics namespace
    
    return point_source_text 

def ArrayOfPointSources(pnt_src_data):
    #create an array of all point source text expressions 
    
    logger.info('Point sources params:')
    logger.info('r')
    logger.info(pnt_src_data.r)
    logger.info('I & disp')
    logger.info(pnt_src_data.i_disp)
    logger.log_n_output("Point sources:", 'magenta')
    
    pnt_src_text = []
    for i in range(len(pnt_src_data.r)):
        pnt_src_text.append(CreatePointSource(pnt_src_data.r[i], pnt_src_data.i_disp[i][0], pnt_src_data.alpha * pnt_src_data.i_disp[i][1]))
        
    return pnt_src_text

def Array_Expression(text_array):
    expression_array = [None]*len(text_array)
    for i in range(len(text_array)):
        expression_array[i] = Expression(text_array[i], degree = 2)
    
    return expression_array

def Contour_plot(r_area, z_area, u, path, f_expr, mesh, plot_title, contour_amount):
    tol, point_num = 0.001, DEFAULT_MESH + 1  # avoid hitting points outside the domain
    r = numpy.linspace(r_area[0] + tol, r_area[1] - tol, point_num)
    z = numpy.linspace(z_area[0] + tol, z_area[1] - tol, int(point_num*mesh[1]/mesh[0]))
    
    levels = numpy.linspace(-15e-5, 0)
    u_contour = numpy.zeros([len(z), len(r)])
    
    for i in range(len(r)):
        for j in range(len(z)):
            u_contour[j, i] = u(r[i], z[j])
    
    if u_contour.max() == u_contour.min():
        print(colored("\u03C8 is the save everywhere!", 'red'))
        print( colored( 'u_max = ', 'green') + str(u.vector()[:].max()) )
        print( colored( 'u_min = ', 'green') + str(u.vector()[:].min()) )
        print( colored( 'u_max - u_min = ', 'green') + str(u.vector()[:].max() - u.vector()[:].min()) )
        
        logger.info( "u_max = %s" % str(u.vector()[:].max()) )
        logger.info( "u_min = %s" % str(u.vector()[:].min()) )
        logger.info( "u_max - u_min = %s" % str(u.vector()[:].max() - u.vector()[:].min()) )
        return 0
    else:
        matplt.contour(r, z, u_contour, levels)
        matplt.xlim(r_area[0], r_area[1])
        matplt.ylim(z_area[0], z_area[1])
        
        matplt.xlabel("r, м")
        matplt.ylabel("z, м")
        matplt.colorbar()
        matplt.gca().set_aspect("equal")
        
    print( colored( 'u_max = ', 'green') + str(u.vector()[:].max()) )
    print( colored( 'u_min = ', 'green') + str(u.vector()[:].min()) )
    print( colored( 'u_max - u_min = ', 'green') + str(u.vector()[:].max() - u.vector()[:].min()) )
    
    logger.info( "u_max = %s" % str(u.vector()[:].max()) )
    logger.info( "u_min = %s" % str(u.vector()[:].min()) )
    logger.info( "u_max - u_min = %s" % str(u.vector()[:].max() - u.vector()[:].min()) )
    
    
    Save_figure(f_expr, mesh[0], mesh[1], '', path, plot_title)
    return 0
            
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
    file = open(file_path, "a") # append write to file mode
    
    text = "%s,%s,%s\n" % (mesh_r, err_max, err_L2)
    file.write(text)
    file.close()
    
    print(colored("Data saved to PATH: %s" % file_path, 'green'))
    
def Plot_error_vs_mesh(name): # u max as a function of mesh parameters on the same solution area
    with open("%svsmesh.txt" % TEXT_FILE_ERROR, "r") as file:
        data = [[float(num) for num in line.split(',')] for line in file]
        
    mesh = Column(data, 0) 
    err_max = Column(data, 1)
    err_L2 = Column(data, 2)
    
    matplt.semilogy(mesh, err_max, 'o', linewidth=2)  # magnify w
    # matplt.legend(["u_max vs default mesh size"], loc='best')
    matplt.grid(True)
    matplt.xlabel('Размер сетки')
    matplt.ylabel('Максимальная ошибка')

    file_path = 'Figures/Post_analyt/analt_errorvsmesh.png'
    matplt.savefig(file_path, dpi = DPI, bbox_inches="tight")
    print(colored("2D plot saved to PATH: %s" % file_path, 'green'))
    
    matplt.close() # close created plot
    
def To_float(arr):
    arr_str = []
    for i in arr:
        arr_str.append(str(i))
        
    return arr_str
    
def Hand_input(p_pow, F_pow):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')
    
    psi = sympy.symbols('u') # flux function #think tomorrow how to define argument psi!
    x = sympy.symbols('x[0]') # r coordinate. used for easy writing of expressions

    p_psi = pow(psi/M.psi_axis, int(p_pow)) #pressure function
    F_psi_2 = 1 - pow(psi/M.psi_axis, int(F_pow)) # poloidal current function

    dp_psi = sympy.diff(p_psi, psi) #pressure and F deriviation
    dF_psi_2 = sympy.diff(F_psi_2, psi) #compiler breaks when 

    f_text = (M0 * pow(x, 2) * M.p_axis * dp_psi + 0.5 * M.F0_2 * dF_psi_2) #right hand expression
    
    f_text = sympy.printing.ccode(f_text)
    f_text = f_text.replace('exp', 'std::exp') # reason being faulty fenics namespace
    # f_text = f_text.replace('pow', 'std::pow') # reason being faulty fenics namespace

    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output(f_text, 'white')

    return f_text

def Initial_guess_for_u(u, const):
    for i in range(len(u.vector())):
        u.vector()[i] = float(const)
    logger.log_n_output_colored_message(colored_message="Initial guess for u: ", color='green', white_message=str(const))
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
        [cell_markers, global_mesh_index] = refine_subdomain(mesh, submesh, cell_markers) #domain numbering starts wwith 1
    
    mesh = refine(mesh, cell_markers)
    logger.info("Mesh refined: %s" % str(global_mesh_index))
    logger.info("After refinement. Number of cells: %d, Number of vertices: %d" % (mesh.num_cells(), mesh.num_vertices()))
    
    return mesh

def refine_subdomain(mesh, submesh, cell_markers):
    bound_box = mesh.bounding_box_tree()
    global_mesh_index = []
    
    for cell in cells(submesh):
        global_mesh_index.append( bound_box.compute_first_entity_collision(cell.midpoint()) )
    
    for cell in cells(mesh):
        if cell.index() in global_mesh_index:
            cell_markers[cell] = True
        else:
            cell_markers[cell] = False
    
    return cell_markers, global_mesh_index

def Create_Subdomain(r, disp, segments):
    circle = mshr.Circle(Point(r[0], r[1]), disp*3, segments=segments)
    logger.info("Created subdomain: %s, r=%e, segments=%d" % (str(r), disp*3, segments))
    
    return circle

def Create_Subdomains(alpha, segments):
    circle_list = []
    ps_d = psd.PointSource(alpha)
    
    for i in range(len(ps_d.r)):
        if abs( ps_d.i_disp[i][0] ) > 1e-6:
            circle_list.append( Create_Subdomain(ps_d.r[i], ps_d.i_disp[i][1] * ps_d.alpha, segments) )
    
    return circle_list

def Set_Subdomains(domain, alpha, segments):
    circle_list = Create_Subdomains(alpha, segments)
    for i in range(len(circle_list)):
        domain.set_subdomain( i+1, circle_list[i] )
        
    return domain, len(circle_list)

def plot_mesh(mesh, path):
    plot(mesh)
    
    matplt.xlabel("r, м")
    matplt.ylabel("z, м")
    matplt.gca().set_aspect("equal")
    
    Save_figure( '', 100, 100, '', path, "" )
    
    return 0


def plasma_sources_coefficients_pow_2(p_correction,F_correction):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')
    
    p_coeff = 2 * M0 * M.p_axis / M.psi_axis**2 * p_correction
    F_2_coeff = -M.F0_2 / M.psi_axis**2 * F_correction

    logger.log_n_output_colored_message(colored_message="p_correction = ", color='green', white_message=str(p_correction))
    logger.log_n_output_colored_message(colored_message="F_2_correction = ", color='green', white_message=str(F_correction))
    
    logger.log_n_output_colored_message(colored_message="p_coeff = ", color='green', white_message=str(p_coeff))
    logger.log_n_output_colored_message(colored_message="F_2_coeff = ", color='green', white_message=str(F_2_coeff))
    
    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output("%s*pow(x[0], 2)*u + %s*u" % (p_coeff, F_2_coeff), 'white')
    
    return p_coeff, F_2_coeff

def plasma_sources_coefficients_pow_2_iteration(p_correction, F_correction, psi_axis):
    M = MEPH.MEPhIST()
    print_colored("MEPhIST data:", 'magenta')
    logger.log_n_output(M.__dict__, 'white')
    
    p_coeff = 2 * M0 * M.p_axis / psi_axis**2 * p_correction
    F_2_coeff = -M.F0_2 / psi_axis**2 * F_correction

    logger.log_n_output_colored_message(colored_message="p_correction = ", color='green', white_message=str(p_correction))
    logger.log_n_output_colored_message(colored_message="F_2_correction = ", color='green', white_message=str(F_correction))
    
    logger.log_n_output_colored_message(colored_message="p_coeff = ", color='green', white_message=str(p_coeff))
    logger.log_n_output_colored_message(colored_message="F_2_coeff = ", color='green', white_message=str(F_2_coeff))
    
    logger.log_n_output("Right hand part: ", 'magenta')
    logger.log_n_output("%s*pow(x[0], 2)*u + %s*u" % (p_coeff, F_2_coeff), 'white')
    
    return p_coeff, F_2_coeff

def spheromak_point(r, R, alpha):
    z = 1/2/alpha * math.sqrt(2 * R**2 - (r - R)**2)
    return z
    
def plot_spheromak_boundary(R, alpha, smoothness):
    r_array = numpy.linspace((1 + math.sqrt(2)) * R, (1 - math.sqrt(2)) * R, smoothness)
    z_array = []
    for r in r_array:
        z_array.append(spheromak_point(r, R, alpha))
    
    z_array = numpy.array(z_array)
    
    r_array = numpy.append(r_array, numpy.flip(r_array)) # замкнуть кривую
    z_array = numpy.append(z_array, numpy.flip(-z_array)) # замкнуть кривую
    return r_array, z_array
    
def spheromak_boundary(R, alpha, smoothness):
    [r_array, z_array] = plot_spheromak_boundary(R, alpha, smoothness)
    boundary_geometry = []
    for i in range(len(r_array)):
        boundary_geometry.append(Point(r_array[i],z_array[i]))
    
    return mshr.Polygon(boundary_geometry)

def spheromak_pressure(psi_0, R, alpha):
    L = Expression("pow(x[0], 2) / pow(%s, 4) * (1 + pow(%s, 2)) * %s" % (R, alpha, psi_0), degree = 2)
    logger.log_n_output("Right hand part:", 'red')
    logger.log_n_output(L._cppcode, 'white')
    return L

def countour_plot_via_mesh(geometry, u, levels, PATH, plot_title):
    u_min = u.vector()[:].min()
    u_max = u.vector()[:].max()
    if u_min == u_max:
        logger.log_n_output(message="Trivial solution. u = %s" % u_max, color='red')
    else:
        triang = tri.Triangulation(*geometry.mesh.coordinates().reshape((-1, 2)).T, triangles=geometry.mesh.cells())
        u_array = u.compute_vertex_values(geometry.mesh)
        
        matplt.tricontour(triang, u_array, levels)
        # matplt.xlim(geometry.r1, geometry.r2)
        # matplt.ylim(geometry.z1, geometry.z2)
        # matplt.xlim(0.05, 0.7)
        # matplt.ylim(-0.6, 0.6)

        matplt.xlabel("r, м")
        matplt.ylabel("z, м")
        matplt.colorbar()
        matplt.gca().set_aspect("equal")
        
        logger.log_n_output_colored_message(colored_message="u_max = ", color='green', white_message=str(u_max))
        logger.log_n_output_colored_message(colored_message="u_min = ", color='green', white_message=str(u_min))
        logger.log_n_output_colored_message(colored_message="u_max-u_min = ", color='green', white_message=str(u_max-u_min))
        
        save_contour_plot(PATH, plot_title)
        
        return u_max
    
def save_contour_plot(PATH, plot_title):
    time_title = Time_name()

    path_my_file = 'Figures/%s/%s' % (PATH, time_title)
    file_path = "%s.png" % path_my_file
    logger.info(file_path)

    matplt.title(plot_title) # titled figure for my self
    matplt.savefig(file_path, dpi = DPI, bbox_inches="tight") #no title figure for reports
    matplt.close() # close created plot
    
    print_colored_n_white(colored_text="3D countour plot saved to PATH: ", color='green', white_text=file_path)
    time.sleep(0.5)

def print_colored_n_white(colored_text, color, white_text):
    print(colored(colored_text, color) + white_text)
    
def ErrorEstimate(u, u_D, mesh):
    error_L2 = errornorm(u_D, u, 'L2')

    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = numpy.max(numpy.abs(vertex_values_u_D - vertex_values_u))

    logger.log_n_output_colored_message(colored_message="error_L2 = ", color='red', white_message=str(error_L2))
    logger.log_n_output_colored_message(colored_message="error_max = ", color='red', white_message=str(error_max))
    
    return error_L2, error_max
# %%

def plot_1D(PATH, u, geometry):
    
    plot(u)
    
    matplt.grid("True")
    matplt.xlim(geometry.a, geometry.b)
    matplt.xlabel("r, м")
    matplt.ylabel("p(r), Па")
    
    logger.print_colored_n_white(colored_text="u_max = ", color='green', white_text=str(u.vector()[:].max()))
    logger.print_colored_n_white(colored_text="u_min = ", color='green', white_text=str(u.vector()[:].min()))
    logger.print_colored_n_white(colored_text="u_max-u_min = ", color='green', white_text=str(u.vector()[:].max() - u.vector()[:].min()))
    
    save_contour_plot(PATH, "")
    
def multiply_u_by_const(u, const):
    for i in range(len(u.vector())):
        u.vector()[i] = float(const)*u.vector()[i]
    logger.log_n_output_colored_message(colored_message="u is multiplied by: ", color='green', white_message=str(const))
    return u
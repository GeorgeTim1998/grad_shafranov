from os import name

from numpy import mat
from imports import *
import time

DPI = 200 # quality of plots
TEXT_FILE_U_MAX = "Text_data/func_max"
TEXT_FILE_ERROR = "Text_data/error"
TEXT_FILE_2D_PLOT = "Text_data/2Dplot"
TWOD_PLOT_SAVE_PATH = 'Figures/Post_analyt'
M0 = 1.25e-6
DEFAULT_MESH = 100

SQ_MIN = 1
SQ_MAX = 2
SQUARE_SIZE_ARRAY = numpy.linspace(SQ_MIN, SQ_MAX, 1+int((SQ_MAX-SQ_MIN)/SQ_MIN))

EPS = 25/13
ALPHA = 0.4
KAPPA = 2

matplt.rcParams["font.family"] = "Times New Roman"

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
    matplt.savefig(file_path, dpi = DPI)
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
    matplt.ylabel('\u03C8}')
    matplt.grid(True)
    matplt.legend()

    file_path = "%s/2D_plots_together_%s.png" % (TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi = 2*DPI)
    matplt.close() # close created plot
    print(colored("2D plots saved together at PATH: %s" % file_path, 'green'))
    
    # now build errors # maybe change it for the better next time!
    iteration = 0
    n = 0
    for i in i_array: 
        matplt.plot(x[n:len(x)-n], delta_arr[iteration][n:len(x)-n], linewidth=0.5, label="Square size = %s" % int(i))
        iteration = iteration + 1
    matplt.xlabel('r')
    matplt.ylabel("\u0394, %")
    matplt.grid(True)
    matplt.legend()
    
    file_path = "%s/2D_plots_together_error_%s.png" % (TWOD_PLOT_SAVE_PATH, int(SQUARE_SIZE_ARRAY[0]))
    matplt.savefig(file_path, dpi = 2*DPI)
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

    # matplt.title("%s - %s\n%s" % (plot_title, mesh_title, f_expr._cppcode)) # titled figure for my self
    
    file_path = "%s%s.png" % (path_my_file, addition)
    matplt.savefig(file_path, dpi = DPI) #no title figure for reports
    matplt.close() # close created plot
    
    print(colored("3D countour plot saved to PATH: %s" % file_path, 'green'))
    
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
    
    matplt.savefig("Figures/umax_vs_mesh_%s.png" % name, dpi = DPI)
    
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
    
    matplt.savefig("Figures/umax_vs_sq_sz_%s_%s.png" % (name, default_mesh_size), dpi = DPI)
    
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

def ErrorEstimate(u, u_D, mesh):
    # Compute error in L2 norm
    error_L2 = errornorm(u_D, u, 'L2')

    # Compute maximum error at vertices
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = numpy.max(numpy.abs(vertex_values_u_D - vertex_values_u))

    # Print errors
    print(colored('error_L2  = ', 'red'), error_L2)
    print(colored('error_max = ', 'red'), error_max)
    
    return error_L2, error_max

def CreatePointSource(r, I, disp):
    x = sympy.symbols('x[0]') # r coordinate
    z = sympy.symbols('x[1]') # r coordinate

    pre_exp = M0/pi/disp/disp * I * x # in sympy write stuff that works
    inner_exp = - (pow(x - r[0], 2) + pow(z - r[1], 2)) / pow(disp, 2) # in sympy write stuff that works
    pre_exp_text = sympy.printing.ccode(pre_exp) # transfer it to text
    inner_exp_text = sympy.printing.ccode(inner_exp) # transfer it to text
    
    point_source_text = "%s*exp(%s)" % (pre_exp_text, inner_exp_text) # assemble function of the point source
    print(colored("Point source: \n", 'magenta') + point_source_text)
    return point_source_text 

def ArrayOfPointSources(pnt_src_data):
    #create an array of all point source text expressions 
    
    pnt_src_text = []
    for i in range(len(pnt_src_data.r)):
        pnt_src_text.append(CreatePointSource(pnt_src_data.r[i], pnt_src_data.i_disp[i][0], pnt_src_data.i_disp[i][1]))
        
    return pnt_src_text

def Array_Expression(text_array):
    expression_array = [None]*len(text_array)
    for i in range(len(text_array)):
        expression_array[i] = Expression(text_array[i], degree = 2)
    
    return expression_array

def Contour_plot(r_area, z_area, u, path, f_expr, mesh, plot_title):
    tol, point_num = 0.001, DEFAULT_MESH + 1  # avoid hitting points outside the domain
    r = numpy.linspace(r_area[0] + tol, r_area[1] - tol, point_num)
    z = numpy.linspace(z_area[0] + tol, z_area[1] - tol, int(point_num*mesh[1]/mesh[0]))
    
    u_contour = numpy.zeros([len(z), len(r)])
    
    for i in range(len(r)):
        for j in range(len(z)):
            u_contour[j, i] = u(r[i], z[j])
            
    matplt.contour(r, z, u_contour)
    matplt.xlim(r_area[0], r_area[1])
    matplt.ylim(z_area[0], z_area[1])
    
    matplt.xlabel("r")
    matplt.ylabel("z")
    matplt.grid()
    matplt.colorbar()
    
    Save_figure(f_expr, mesh[0], mesh[1], 'cont_title', path, plot_title)
            
def Mesh_to_xml():
    file = File('Mesh/file.xml')

def D_config(smoothness):
    tol = 0.0001
    t = numpy.linspace(tol, math.pi - tol, smoothness)
    x = numpy.ones(smoothness) + EPS*numpy.cos(t + ALPHA*numpy.sin(t))
    z = EPS*KAPPA*numpy.sin(t)
    
    x = numpy.append(x, numpy.flip(x)) # move it along x axis!
    z = numpy.append(z, numpy.flip(-z))
    matplt.plot(x, z)
    matplt.grid()
    matplt.show()
    
def Write2file_errors(mesh_r, mesh_z, err_L2, err_max):
    file_path = "%s_%s.txt" % (TEXT_FILE_ERROR, mesh_r)
    file = open(file_path, "a") # append write to file mode
    
    text = "%s,%s,%s\n" % (mesh_r, err_max, err_L2)
    file.write(text)
    file.close()
    
    print(colored("Data saved to PATH: %s" % file_path, 'green'))
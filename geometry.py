from fenics import *
import logger
import mshr
import numpy
class Geometry:
    def operator_weights(self, V):
        r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
        r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such
        
        return r_2, r
#%% Rectangle geometry
    def rectangle_mesh_init(self, r1, r2, z1, z2, default_mesh):
        self.r1 = r1
        self.r2 = r2
        self.z1 = z1
        self.z2 = z2
        self.default_mesh = default_mesh
        
        self.rectangle_mesh()
        self.rectangle_mesh_log()
        
    def rectangle_area(self):
        self.area = [self.r1, self.r2, self.z1, self.z2]
        rect_low = Point(self.area[0], self.area[2]) #define rectangle size: lower point
        rect_high = Point(self.area[1], self.area[3]) #define rectangle size: upper point
        
        return rect_low, rect_high
        
    def rectangle_mesh(self):
        [rect_low, rect_high] = self.rectangle_area()
        self.rectangle_mesh_values()
        
        self.mesh = RectangleMesh(rect_low, rect_high, self.mesh_r, self.mesh_z) # points define domain size rect_low x rect_high
        
    def rectangle_mesh_log(self):
        logger.log_n_output_colored_message(colored_message="mesh_r = ", color='green', white_message=str(self.mesh_r))
        logger.log_n_output_colored_message(colored_message="mesh_z = ", color='green', white_message=str(self.mesh_z))
        logger.info('R1 = %f, Z1 = %f' % (self.r1, self.z1))
        logger.info('R2 = %f, Z2 = %f' % (self.r2, self.z2))
        logger.info( "Number of cells: %d, Number of vertices: %d" % (self.mesh.num_cells(), self.mesh.num_vertices()) )
    
    def rectangle_mesh_values(self):
        self.mesh_r, self.mesh_z = self.default_mesh, abs(int(self.default_mesh * (self.z2-self.z1)/(self.r2-self.r1)))
        
#%% Arbitrary mesh
    def arbitrary_mesh_init(self):
        domain = mshr.Circle(Point(0.3, 0), 0.25)
        self.mesh = mshr.generate_mesh(domain, 50)
#%% 1D geometry (Interval geometry)
    def interval_mesh_init(self, a, b, default_mesh):
        self.a = a
        self.b = b
        self.default_mesh = default_mesh
        
        self.mesh = IntervalMesh(default_mesh, a,b)
        self.interval_mesh_log()
        
    def interval_mesh_log(self):
        logger.log_n_output_colored_message(colored_message="a = ", color='green', white_message=str(self.a))
        logger.log_n_output_colored_message(colored_message="b = ", color='green', white_message=str(self.b))
        logger.log_n_output_colored_message(colored_message="default mesh = ", color='green', white_message=str(self.default_mesh))

#%% Create domains
    def rectangle_domain(self, area):
        self.r1 = area[0]
        self.r2 = area[1]
        self.z1 = area[2]
        self.z2 = area[3]
        
        self.domain_centre_r = (area[0] + area[1]) / 2
        self.domain_centre_z = (area[2] + area[3]) / 2
        
        rect_low = Point(self.r1, self.z1) #define rectangle size: lower point
        rect_high = Point(self.r2, self.z2) #define rectangle size: upper point
        
        self.rect_domain = mshr.Rectangle(rect_low, rect_high)
        
        self.rectangle_domain_log()

        return self.rect_domain
    
    def rectangle_domain_log(self):
        logger.info('R1 = %f, Z1 = %f' % (self.r1, self.z1))
        logger.info('R2 = %f, Z2 = %f' % (self.r2, self.z2))
        
    def circle_domain(self, radius, segments):
        circle = mshr.Circle(Point(self.domain_centre_r, self.domain_centre_z), radius, segments=segments)
        logger.info("Created subdomain: centre: %s, radius=%f, segments=%d" % (str([self.domain_centre_r, self.domain_centre_z]), radius, segments))
        
        return circle
    
    def generate_mesh_in_domain(self, domain, density):
        self.density = density
        self.mesh = mshr.generate_mesh(domain, density)
        self.log_mesh_in_domain()
        
    def log_mesh_in_domain(self):
        logger.log_n_output_colored_message(colored_message="Mesh density = ", color='green', white_message=str(self.density))
        logger.info( "Number of cells: %d, Number of vertices: %d" % (self.mesh.num_cells(), self.mesh.num_vertices()) )
#%% MEPhIST domain
    def mephist_vessel(self):
        [x, z] = self.read_data_from_file()
        
        point_list = []

        for i in range(len(x)):
            point_list.append(Point(x[i], z[i]))
            
        self.domain = mshr.Polygon(point_list)
        
        return self.domain
    
    def Column(self, matrix, col):
        return [row[col] for row in matrix]
    
    def read_data_from_file(self):
        file_path = "КонтурКамерыСредний.txt"
        with open(file_path, "r") as file: # change to Read_from_file func
            data = [[float(num) for num in line.split('  ')] for line in file]

        x = numpy.array(self.Column(data, 0))*1e-3
        z = numpy.array(self.Column(data, 1))*1e-3
        
        camera_height = 581.69*1e-3
        z_height = z[-2]

        x = x[0:len(x)-1]
        z = z[0:len(z)-1]

        z = z - z_height/2 * numpy.ones(len(z))
        z = z + (camera_height/2 - numpy.amax(z)) * numpy.ones(len(z)) 
        
        x = numpy.append(x, numpy.flip(x)) # move it along x axis!
        z = numpy.append(z, numpy.flip(-z))
        
        return x, z
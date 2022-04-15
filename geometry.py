from fenics import *
import logger

class Geometry:
    def rectangle_mesh_init(self, r1, r2, z1, z2, default_mesh):
        self.r1 = r1
        self.r2 = r2
        self.z1 = z1
        self.z2 = z2
        self.default_mesh = default_mesh
        
    def rectangle_area(self):
        self.area = [self.r1, self.r2, self.z1, self.z2]
        rect_low = Point(self.area[0], self.area[2]) #define rectangle size: lower point
        rect_high = Point(self.area[1], self.area[3]) #define rectangle size: upper point
        
        return rect_low, rect_high
    
    def rectangle_mesh_values(self):
        self.mesh_r, self.mesh_z = self.default_mesh, abs(int(self.default_mesh * (self.z2-self.z1)/(self.r2-self.r1)))
        
    def rectangle_mesh(self):
        [rect_low, rect_high] = self.rectangle_area(self)
        self.rectangle_mesh_values()
        
        self.mesh = RectangleMesh(rect_low, rect_high, self.mesh_r, self.mesh_z) # points define domain size rect_low x rect_high
        
        return self.mesh
        
    def rectangle_mesh_log(self):
        logger.log_n_output('DEFAULT_MESH = %d' % self.default_mesh, 'green')
        logger.info('R1 = %f, Z1 = %f' % (self.r1, self.z1))
        logger.info('R2 = %f, Z2 = %f' % (self.r2, self.z2))
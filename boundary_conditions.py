from fenics import *
import sympy
import logger

class BoundaryConditions:
#%% Spheromak
    def spheromak_boundary_condition(self, psi_0, R, alpha):
        self.psi_0 = psi_0
        self.R = R
        self.alpha = alpha
        
        self.spheromak_sympy()
        self.spheromak_source()
        self.log_spheromak_boundary()
        self.psi_sol_expr = Expression(self.spheromak_text, degree = 2)

    def spheromak_sympy(self):
        r = sympy.symbols('x[0]')
        z = sympy.symbols('x[1]')
        
        symbols = self.psi_0 * pow(r, 2) / pow(self.R, 4) * ( 2*pow(self.R, 2) - pow(r, 2) - 4*pow(self.alpha, 2)*pow(z, 2) )
        self.spheromak_text = sympy.printing.ccode(symbols)
    
    def spheromak_source(self):
        r = sympy.symbols('x[0]')
        symbols = 8 * self.psi_0 * pow(r, 2) / pow(self.R, 4) * (1 + pow(self.alpha, 2))
        self.spheromak_right_hand_expr = Expression(sympy.printing.ccode(symbols), degree = 2)
        
        return self.spheromak_right_hand_expr
        
    def log_spheromak_boundary(self):
        logger.info("R = %s, alpha = %s" % (self.R, self.alpha))
        logger.info("psi_sol = %s" % self.spheromak_text)
        logger.log_n_output("Rigth hand part:", 'red')
        logger.log_n_output(self.spheromak_right_hand_expr._cppcode, 'white')
#%% Constant boundary conditions
    def constant_boundary_condition(self, const_str):
        logger.info("u_D = %s" % const_str)
        self.const_boundary = Expression(const_str, degree = 1) # Define boundary condition
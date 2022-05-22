# Файл, в котором хранится информация об граничных условиях для задач, допускающих аналитическое решение
# А также граничные условия с постоянным значением функции на границе также задаются здесь

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
        logger.info("psi_0 = %s" % self.psi_0)
        logger.info("psi_sol = %s" % self.spheromak_text)
        logger.log_n_output("Rigth hand part:", 'red')
        logger.log_n_output(self.spheromak_right_hand_expr._cppcode, 'white')
#%% Constant boundary conditions
    def constant_boundary_condition(self, const_str):
        logger.info("u_D = '%s'" % const_str)
        self.const_boundary_condition = Expression(const_str, degree = 1) # Define boundary condition
        
        return self.const_boundary_condition
#%% Tokamak D configuration
    def tokamak_D_config_boundary_condition(self, psi_0, eps, a, b):
        self.psi_0 = psi_0
        self.eps = eps
        self.a = a
        self.b = b
        
        self.tokamak_D_config_sympy()
        self.tokamak_D_config_source()
        self.log_tokamak_D_config_boundary()
        self.psi_sol_expr = Expression(self.tokamak_D_config_text, degree = 2)

    def tokamak_D_config_sympy(self):
        r = sympy.symbols('x[0]')
        z = sympy.symbols('x[1]')
        
        symbols = self.psi_0 * (pow(r, 2) / pow(self.eps, 2) - 1) * (1 - pow(r, 2) / pow(self.a, 2) - pow(z, 2) / pow(self.b, 2))
        self.tokamak_D_config_text = sympy.printing.ccode(symbols)
    
    def tokamak_D_config_source(self):
        r = sympy.symbols('x[0]')
        symbols = 2*self.psi_0*( (pow(2/self.a, 2) + pow(1/self.b, 2)) * pow(r/self.eps, 2) - pow(1/self.b, 2) )
        self.tokamak_D_config_right_hand_expr = Expression(sympy.printing.ccode(symbols), degree = 2)
        
        return self.tokamak_D_config_right_hand_expr
        
    def log_tokamak_D_config_boundary(self):
        logger.log_n_output_colored_message(colored_message="psi_0 = ", color = 'green', white_message=str(self.psi_0))
        
        logger.log_n_output_colored_message(colored_message="eps (ellipse) = ", color = 'green', white_message=str(self.eps))
        logger.log_n_output_colored_message(colored_message="a = ", color = 'green', white_message=str(self.a))
        logger.log_n_output_colored_message(colored_message="b = ", color = 'green', white_message=str(self.b))
        
        logger.info("psi_sol = %s" % self.tokamak_D_config_text)
        logger.log_n_output("Rigth hand part:", 'red')
        logger.log_n_output(self.tokamak_D_config_right_hand_expr._cppcode, 'white')
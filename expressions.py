from fenics import *
import sympy
import logger

M0 = 1.25e-6

class Expressions:
#%% Axisymmetric configuration
    def axissymm_config(self, j0, r0):
        self.j0 = j0
        self.r0 = r0
        
        r = sympy.symbols('x[0]')
        
        f_symbols = -M0 * pow(j0, 2) * pow(r0/r, 2)
        self.axissymm_config_right_hand_expr = Expression(sympy.printing.ccode(f_symbols), degree = 2)
        self.axissymm_config_log()
        self.axissymm_config_solution()
        
    def axissymm_config_log(self):
        logger.log_n_output_colored_message(colored_message='j0 = ', color='green', white_message=str(self.j0))
        logger.log_n_output_colored_message(colored_message='r0 = ', color='green', white_message=str(self.r0))
        
        logger.print_colored("Right hand part: ", 'magenta')
        logger.print_colored(self.axissymm_config_right_hand_expr._cppcode, 'white')
    
    def axissymm_config_solution(self):
        r = sympy.symbols('x[0]')
        
        p_symbols = M0 * pow(self.j0, 2) * pow(self.r0, 2) * (1/r - 1/self.r0)
        self.axissymm_config_solution_expr = Expression(sympy.printing.ccode(p_symbols), degree = 2)
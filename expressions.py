from fenics import *
import sympy
import logger
import funcs as fu
from MEPHIST_data import MEPhIST
import math

M0 = 1.25e-6


class Expressions:
    # %% Axisymmetric configuration
    def axissymm_config(self, j0, r0):
        self.j0 = j0
        self.r0 = r0

        r = sympy.symbols('x[0]')

        f_symbols = -M0 * pow(j0, 2) * pow(r0/r, 2)
        self.axissymm_config_right_hand_expr = Expression(
            sympy.printing.ccode(f_symbols), degree=2)
        self.axissymm_config_log()
        self.axissymm_config_solution()

    def axissymm_config_log(self):
        logger.log_n_output_colored_message(
            colored_message='j0 = ', color='green', white_message=str(self.j0))
        logger.log_n_output_colored_message(
            colored_message='r0 = ', color='green', white_message=str(self.r0))

        logger.print_colored("Right hand part: ", 'magenta')
        logger.print_colored(
            self.axissymm_config_right_hand_expr._cppcode, 'white')

    def axissymm_config_solution(self):
        r = sympy.symbols('x[0]')

        p_symbols = M0 * pow(self.j0, 2) * pow(self.r0, 2) * (1/r - 1/self.r0)
        self.axissymm_config_solution_expr = Expression(
            sympy.printing.ccode(p_symbols), degree=2)

# %% Move  Spheromak source

    def moving_sphmk_source(self, R, a, alpha, psi0, t, problem):
        logger.log_n_output_colored_message(
            colored_message="Displacement: ",
            color='green',
            white_message="disp_fact*vessl_rad*t/tm = %s*%s*%.3e/%.3e" % (str(problem.disp_fact), str(problem.ves_inner_radius), t, problem.tm))
        fu.print_colored_n_white(colored_text='ts = ',
                                 color='green', white_text="%.3e" % problem.ts)
        logger.log_n_output_colored_message(
            colored_message="Displacement of source by: ", color='green', white_message=str(a))
        x = sympy.symbols('x[0]')
        z = sympy.symbols('x[1]')

        expr = psi0 * (x-a)**2 * (2*R**2 - (x-a)**2 -
                                  4*alpha**2*(z-0.1)**2) / R**4
        # expr = psi0*(R**2 - (x-R-a)**2 - (z-0.3)**2) / R**2

        # first term (смотри формы оператора ГШ)
        f_expr_x1 = sympy.diff(sympy.diff(expr, x), x)
        f_expr_x2 = -1/x*sympy.diff(expr, x)
        f_expr_z = sympy.diff(sympy.diff(expr, z), z)

        f_expr = -sympy.simplify(f_expr_x1 + f_expr_x2 + f_expr_z)

        f_text = sympy.printing.ccode(f_expr)

        fu.print_colored(text='Moving source:', color='green')
        sympy.pprint(f_expr)

        return Expression(f_text, degree=2)

    def moving_point_source(self, R, a, t, problem):
        M = MEPhIST()

        logger.log_n_output_colored_message(
            colored_message="Displacement: ",
            color='green',
            white_message="disp_fact*vessl_rad*t/tm = %s*%s*%.3e/%.3e" % (str(problem.disp_fact), str(problem.ves_inner_radius), t, problem.tm))
        fu.print_colored_n_white(colored_text='ts = ',
                                 color='green', white_text="%.3e" % problem.ts)
        logger.log_n_output_colored_message(
            colored_message="Displacement of source by: ", color='green', white_message=str(a))

        x = sympy.symbols('x[0]')
        z = sympy.symbols('x[1]')

        sigma = problem.ves_inner_radius*0.25
        j0 = M.I/(math.pi*sigma**2)
        logger.info(message="sigma = %.3e" % sigma)
        logger.info(message="j0 = I/(pi*sigma**2) = %.3e" % j0)

        f_expr = M0*x*j0 * sympy.exp(-((x-R-a)**2 + z**2) / sigma**2)

        f_text = sympy.printing.ccode(f_expr)
        f_text = f_text.replace('exp', 'std::exp')

        fu.print_colored(text='Moving source:', color='green')
        sympy.pprint(f_expr)

        return Expression(f_text, degree=2)
    
    def point_source_t0(self, R, problem):
        M = MEPhIST()

        x = sympy.symbols('x[0]')
        z = sympy.symbols('x[1]')

        sigma = problem.ves_inner_radius*0.25
        j0 = M.I/(math.pi*sigma**2)
        logger.info(message="sigma = %.3e" % sigma)
        logger.info(message="j0 = I/(pi*sigma**2) = %.3e" % j0)

        f_expr = M0*x*j0 * sympy.exp(-((x-R)**2 + z**2) / sigma**2)

        f_text = sympy.printing.ccode(f_expr)
        f_text = f_text.replace('exp', 'std::exp')

        fu.print_colored(text='Initial source:', color='green')
        sympy.pprint(f_expr)

        return Expression(f_text, degree=2)

import sympy as s
import IPython.display as disp

x = s.symbols('x[0]')
z = s.symbols('x[1]')

r0 = s.symbols('r0')
R = s.symbols('R')
a = s.symbols('a')

expr = (x-r0)**2 * (2*R**2 - (x-r0)**2 - 4*a**2*z**2) / R**4

f_text = s.printing.ccode(expr)
# print(f_text)

f_expr_x1 = s.diff(s.diff(expr, x), x) # first term (смотри формы оператора ГШ)
f_expr_x2 = -1/x*s.diff(expr, x) 
f_expr_z = s.diff(s.diff(expr, z), z)

f_expr = -s.simplify(f_expr_x1 + f_expr_x2 + f_expr_z)

s.pprint(f_expr)
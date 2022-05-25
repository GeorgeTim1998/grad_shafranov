import sympy as s

x = s.symbols('x')  # r coordinatez
z = s.symbols('z')
R = s.symbols('R')
a = s.symbols('a')

expr = x**2 * (2*R**2 - x**2 - 4*a**2*z**2) / R**4

f_text = s.printing.ccode(expr)
print(f_text)

f_expr_x1 = s.diff(s.diff(expr, x), x) # first term (смотри формы оператора ГШ)
f_expr_x2 = -1/x*s.diff(expr, x)

f_expr_z = s.diff(s.diff(expr, z), z)

f_expr = s.simplify(f_expr_x1 + f_expr_x2 + f_expr_z)
print(s.printing.ccode(f_expr))
# print(s.printing.ccode(f_text_x))
# print(s.printing.ccode(f_text_z))
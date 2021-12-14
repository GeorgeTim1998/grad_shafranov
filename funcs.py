from __future__ import print_function
from math import degrees
from fenics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, gcf, isinteractive, title
import datetime
import numpy 
import sympy
from termcolor import colored
import pylab as plt

def Form_f_text(A1, A2):
    #A1 = 4*mo*p', A2 = FF'
    #deriviation are calculated using sympy library
    x = sympy.symbols('x[0]') # r coordinate
    f_text = sympy.printing.ccode(A1 * pow(x, 2) + A2)
    print(colored("INVERCED right-hand equation side: \n", 'magenta') + f_text)

    return f_text
import time

from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy.utilities.lambdify import lambdify
from sympy import diff, Symbol, sin, cos, pi, Matrix, sqrt, acos
from sympy import Eq, solve, trigsimp, simplify
from sympy import init_printing, pprint, pretty
import numpy as np
import dill
from random import random
import pdb

from constants import *

DISCRIMINANT_EPSILON = 1e-10

def logDuration(duration_name, tic):
    if __name__ == "__main__":
        print(str(duration_name) + ": " + str(time.time() - tic) + "s")

def findTheta1(theta2, theta3, x):
    s2 = sin(theta2)
    s3 = sin(theta3)
    s23 = sin(theta2+theta3)
    c2 = cos(theta2)
    c3 = cos(theta3)
    c23 = cos(theta2+theta3)
    # Given theta2, theta3 and x, find theta1
    # i.e. solve                                          l_1*c1 + l_2*c12 + l_3*c123 = x
    # Use sum of angles formula: l_1*c1 + l_2*(c1*c2 - s1*s2) + l_3*(c1*c23 - s1*s23) = x
    # Rearrange terms:                                                q_s*s1 + q_c*c1 = x
    q_s = -l_2*s2 - l_3*s23
    q_c = l_1 + l_2*c2 + l_3*c23
    # Rearrange terms:                                                 q_s*s1 = x - q_c*c1
    # Square both sides:                                         q_s**2*s1**2 = x**2 - 2*x*q_c*c1 + q_c**2*c1**2
    # Rearrange terms    (q_c**2 + q_s**2)*c1**2 - 2*x*q_c*c1 + x**2 - q_s**2 = 0
    # Solve quadratic equation for c1
    a = q_c**2 + q_s**2
    b = -2*x*q_c
    c = x**2 - q_s**2

    discriminant = b**2 - 4*a*c
    # To deal with numerical issues
    # if abs(discriminant) < DISCRIMINANT_EPSILON:
        # discriminant = 0.0
    c1_1 = (-b + sqrt(discriminant))/(2*a)
    c1_2 = (-b - sqrt(discriminant))/(2*a)
    theta1_1 = acos(c1_1)
    theta1_2 = acos(c1_2)
    # return theta1_1, theta1_2
    # theta1 should be in the range of 0 to -pi
    return -theta1_1
    # return theta1_2

def findFrontWheelPosition(theta1, theta2, theta3, is_symbolic = True):
    if is_symbolic:
        sin = pydrake.symbolic.sin
        cos = pydrake.symbolic.cos
    else:
        sin = np.sin
        cos = np.cos

    theta12 = theta2 + theta1
    theta123 = theta3 + theta12
    s1 = sin(theta1)
    s2 = sin(theta2)
    s3 = sin(theta3)
    s12 = sin(theta12)
    s123 = sin(theta123)
    c1 = cos(theta1)
    c2 = cos(theta2)
    c3 = cos(theta3)
    c12 = cos(theta12)
    c123 = cos(theta123)

    x = l_1*c1 + l_2*c12 + l_3*c123
    y = -l_1*s1 - l_2*s12 - l_3*s123
    return (x, y)

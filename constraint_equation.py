from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy.utilities.lambdify import lambdify
from sympy import diff, Symbol, MatrixSymbol, sin, cos, pi, Matrix, sqrt, acos
from sympy import Eq, solve, trigsimp, simplify
import numpy as np
import pdb

from constants import *

# symbolic_trig = [{
        # 'sin': pydrake.symbolic.sin,
        # 'cos': pydrake.symbolic.cos}, 'numpy']
numpy_trig = ['numpy']

t = Symbol('t')

theta1 = dynamicsymbols('theta1')
theta1_d = diff(theta1, t)
theta1_dd = diff(theta1_d, t)
theta2 = dynamicsymbols('theta2')
theta2_d = diff(theta2, t)
theta2_dd = diff(theta2_d, t)
theta3 = dynamicsymbols('theta3')
theta3_d = diff(theta3, t)
theta3_dd = diff(theta3_d, t)
theta4 = dynamicsymbols('theta4')
theta4_d = diff(theta4, t)
theta4_dd = diff(theta4_d, t)

phi = dynamicsymbols('phi')
phi_d = diff(phi, t)
phi_dd = diff(phi_d, t)

q = Matrix([theta1, theta2, theta3, theta4, phi])
q_d = Matrix([theta1_d, theta2_d, theta3_d, theta4_d, phi_d])

theta1_0 = Symbol('theta1_0')
theta1_d_0 = Symbol('theta1_d_0')
theta1_dd_0 = Symbol('theta1_dd_0')
theta2_0 = Symbol('theta2_0')
theta2_d_0 = Symbol('theta2_d_0')
theta2_dd_0 = Symbol('theta2_dd_0')
theta3_0 = Symbol('theta3_0')
theta3_d_0 = Symbol('theta3_d_0')
theta3_dd_0 = Symbol('theta3_dd_0')
theta4_0 = Symbol('theta4_0')
theta4_d_0 = Symbol('theta4_d_0')
theta4_dd_0 = Symbol('theta4_dd_0')
phi_0 = Symbol('phi_0')
phi_d_0 = Symbol('phi_d_0')
phi_dd_0 = Symbol('phi_dd_0')

# ORDER IS IMPORTANT!
# We want to substitute theta2_d before theta2 because Derivative(theta2, t) != Derivative(theta2_0, t)
# Substitution will fail if we substituted theta2 first
substitutions = [
    (theta1_dd, theta1_dd_0),
    (theta2_dd, theta2_dd_0),
    (theta3_dd, theta3_dd_0),
    (theta4_dd, theta4_dd_0),
    (phi_dd, phi_dd_0),
    (theta1_d, theta1_d_0),
    (theta2_d, theta2_d_0),
    (theta3_d, theta3_d_0),
    (theta4_d, theta4_d_0),
    (phi_d, phi_d_0),
    (theta1, theta1_0),
    (theta2, theta2_0),
    (theta3, theta3_0),
    (theta4, theta4_0),
    (phi, phi_0)]

theta_dd_0 = [theta1_dd_0, theta2_dd_0, theta3_dd_0]

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

h = Matrix([
    (findTheta1(theta2, theta3, STEP_POSITION - w_r) - theta1)**2,
    0.0,
    0.0,
    0.0,
    0.0])

H = h.jacobian(q)
calc_H = lambdify(
        [q],
        H)

H_d = H.diff(t)
calc_H_d = lambdify(
        [q, q_d],
        H_d)

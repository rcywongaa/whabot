import time

tic = time.time()
from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy.utilities.lambdify import lambdify
from sympy import diff, Symbol, sin, cos, pi, Matrix
from sympy import Eq, solve, trigsimp
from sympy import init_printing, pprint, pretty
import_sympy_duration = time.time() - tic
import numpy as np

from constants import *
import pdb

import pydrake.symbolic

symbolic_trig = [{
        'sin': pydrake.symbolic.sin,
        'cos': pydrake.symbolic.cos}, 'numpy']
numpy_trig = ['numpy']

with open("res/constants.yaml", 'r') as stream:
    constants = yaml.safe_load(stream)

g = 9.81

tic = time.time()
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
theta = Matrix([theta1, theta2, theta3])
theta_d = Matrix([theta1_d, theta2_d, theta3_d])
theta_dd = Matrix([theta1_dd, theta2_dd, theta3_dd])

theta1_0 = Symbol('theta1_0')
theta1_d_0 = Symbol('theta1_d_0')
theta1_dd_0 = Symbol('theta1_dd_0')
theta2_0 = Symbol('theta2_0')
theta2_d_0 = Symbol('theta2_d_0')
theta2_dd_0 = Symbol('theta2_dd_0')
theta3_0 = Symbol('theta3_0')
theta3_d_0 = Symbol('theta3_d_0')
theta3_dd_0 = Symbol('theta3_dd_0')
theta_0 = [theta1_0, theta2_0, theta3_0]
theta_d_0 = [theta1_d_0, theta2_d_0, theta3_d_0]

# ORDER IS IMPORTANT!
# We want to substitute theta1_d before theta1 because Derivative(theta1, t) != Derivative(theta1_0, t)
# Substitution will fail if we substituted theta1 first
substitutions = [
    (theta1_dd, theta1_dd_0),
    (theta2_dd, theta2_dd_0),
    (theta3_dd, theta3_dd_0),
    (theta1_d, theta1_d_0),
    (theta2_d, theta2_d_0),
    (theta3_d, theta3_d_0),
    (theta1, theta1_0),
    (theta2, theta2_0),
    (theta3, theta3_0)]

theta_dd_0 = [theta1_dd_0, theta2_dd_0, theta3_dd_0]

F = Symbol('F')
tau1 = Symbol('tau1')
tau2 = Symbol('tau2')
tau3 = Symbol('tau3')
# tau4 = Symbol('tau4')

N = ReferenceFrame('N')
i = N.x
j = N.y
k = N.z

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

P2 = l_1*c1*i - l_1*s1*k
P3 = P2 + l_2*c12*i - l_2*s12*k
P4 = P3 + l_3*c123*i - l_3*s123*k
Pb = P2 + l_2/2.0*c12*i - l_2/2.0*s12*k + l_b*cos(theta1 + theta2 + pi/2.0)*i - l_b*sin(theta1 + theta2 + pi/2.0)*k

v2 = P2.diff(t, N)
v3 = P3.diff(t, N)
v4 = P4.diff(t, N)
vb = Pb.diff(t, N)

KE = 0.5*m_2*v2.dot(v2) + 0.5*m_3*v3.dot(v3) + 0.5*m_4*v4.dot(v4) + 0.5*m_b*vb.dot(vb)
PE = m_2*g*P2.dot(k) + m_3*g*P3.dot(k) + m_4*g*P4.dot(k) + m_b*g*Pb.dot(k)
L = Matrix([KE - PE])

lhs = (L.jacobian(theta_d).diff(t) - L.jacobian(theta)).T


# External forces
# FIXME
# F_x = -F_t[0] # Normal reaction of wall
F_x = 0.0*i
F_z = F*k # Vertical force input force (from wheel)
force = F_x + F_z

rhs1 = tau1 + force.dot(P4.diff(theta1, N))
rhs2 = tau2 + force.dot(P4.diff(theta2, N))
rhs3 = tau3 + force.dot(P4.diff(theta3, N))
rhs = Matrix([rhs1, rhs2, rhs3])

eom = Eq(lhs, rhs).subs(substitutions)
formulate_lagrange_duration = time.time() - tic

# tic = time.time()
# eom = trigsimp(eom, method="fu") # Simplify by minimizing trigonometric functions
# simplify_duration = time.time() - tic

tic = time.time()
theta_dd_eom = solve(eom, theta_dd_0, dict=True, simplify=False, rational=False)
solve_eom_duration = time.time() - tic

tic = time.time()
theta1_dd_eom = theta_dd_eom[0][theta1_dd_0]
theta2_dd_eom = theta_dd_eom[0][theta2_dd_0]
theta3_dd_eom = theta_dd_eom[0][theta3_dd_0]

lambda_parameters = [theta_0, theta_d_0, (tau1, tau2, tau3, F)]
theta1_dd_np = lambdify(
        lambda_parameters,
        theta1_dd_eom,
        modules=numpy_trig)
theta1_dd_sym = lambdify(
        lambda_parameters,
        theta1_dd_eom,
        modules=symbolic_trig)

theta2_dd_np = lambdify(
        lambda_parameters,
        theta2_dd_eom,
        modules=numpy_trig)
theta2_dd_sym = lambdify(
        lambda_parameters,
        theta2_dd_eom,
        modules=symbolic_trig)

theta3_dd_np = lambdify(
        lambda_parameters,
        theta3_dd_eom,
        modules=numpy_trig)
theta3_dd_sym = lambdify(
        lambda_parameters,
        theta3_dd_eom,
        modules=symbolic_trig)

lambdify_duration = time.time() - tic

def calc_theta1_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F, is_symbolic = False):
    if is_symbolic:
        theta1_dd_lambd = theta1_dd_sym
    else:
        theta1_dd_lambd = theta1_dd_np
    return theta1_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

def calc_theta2_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F, is_symbolic = False):
    if is_symbolic:
        theta2_dd_lambd = theta2_dd_sym
    else:
        theta2_dd_lambd = theta2_dd_np
    return theta2_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

def calc_theta3_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F, is_symbolic = False):
    if is_symbolic:
        theta3_dd_lambd = theta3_dd_sym
    else:
        theta3_dd_lambd = theta3_dd_np
    return theta3_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

def findTheta1(theta2, theta3, theta4, is_symbolic):
    if is_symbolic:
        sin = pydrake.symbolic.sin
        cos = pydrake.symbolic.cos
    else:
        sin = np.sin
        cos = np.cos

    y = theta4*w_r
    s2 = sin(theta2)
    s3 = sin(theta3)
    s23 = sin(theta2+theta3)
    c2 = cos(theta2)
    c3 = cos(theta3)
    c23 = cos(theta2+theta3)
    # Given theta2, theta3 and y, find theta1
    # i.e. solve                                          l_1*s1 + l_2*s12 + l_3*s123 = y
    # Use sum of angles formula: l_1*s1 + l_2*(s1*c2 + c1*s2) + l_3*(s1*c23 + c1*s23) = y
    # Rearrange terms:                                                q_s*s1 + q_c*c1 = y
    q_s = l_1 + l_2*c2 + l_3*c23
    q_c = -l_2*s2 - l_3*s23
    # Rearrange terms:                                                 q_s*s1 = y - q_c*c1
    # Square both sides:                                         q_s**2*s1**2 = y**2 - 2*y*q_c*c1 + q_c**2*c1**2
    # Rearrange terms    (q_c**2 + q_s**2)*c1**2 - 2*y*q_c*c1 + y**2 - q_s**2 = 0
    # Solve quadratic equation for c1
    a = q_c**2 + q_s**2
    b = -2*y*q_c
    c = y**2 - q_s**2

    discriminant = b**2 - 4*a*c
    # To deal with numerical issues
    if abs(discriminant) < DISCRIMINANT_EPSILON:
        discriminant = 0.0
    c1_1 = (-b + np.sqrt(discriminant))/(2*a)
    c1_2 = (-b - np.sqrt(discriminant))/(2*a)
    theta1_1 = np.arccos(c1_1)
    theta1_2 = np.arccos(c1_2)
    return theta1_1, theta1_2

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

if __name__ == "__main__":
    print("Imported sympy in " + str(import_sympy_duration) + "s")
    print("Lagrange formulated in " + str(formulate_lagrange_duration) + "s")
    # print("Simplify in " + str(simplify_duration) + "s")
    print("EOM solved in " + str(solve_eom_duration) + "s")
    print("Lambdified in " + str(lambdify_duration) + "s")

    from random import random

    tic = time.time()
    for i in range(100):
        calc_theta1_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random(), random())
        calc_theta2_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random(), random())
        calc_theta3_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random(), random())
    print("300 EOM calculation in " + str(time.time() - tic) + "s")

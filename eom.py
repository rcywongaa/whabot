import time

from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy.utilities.lambdify import lambdify
from sympy import diff, Symbol, sin, cos, pi, Matrix, sqrt, acos
from sympy import Eq, solve, trigsimp, simplify, cse
from sympy import init_printing, pprint, pretty
import numpy as np
import pdb

from constants import *
# g = Symbol('g')
# m_w = Symbol('m_w')
# m_1 = Symbol('m_1')
# m_2 = Symbol('m_2')
# m_3 = Symbol('m_3')
# m_4 = Symbol('m_4')
# m_b = Symbol('m_b')
# l_1 = Symbol('l_1')
# l_2 = Symbol('l_2')
# l_3 = Symbol('l_3')
# l_b = Symbol('l_b')
# w_r = Symbol('w_r')
# STEP_DEPTH = Symbol('STEP_DEPTH')
# STEP_WIDTH = Symbol('STEP_WIDTH')
# STEP_HEIGHT = Symbol('STEP_HEIGHT')
# STEP_POSITION = Symbol('STEP_POSITION')

import pydrake.symbolic

symbolic_trig = [{
        'sin': pydrake.symbolic.sin,
        'cos': pydrake.symbolic.cos}, 'numpy']
numpy_trig = ['numpy']

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
    # Rearrange terms:                                                 q_s*c1 = x - q_c*c1
    # Square both sides:                                         q_s**2*c1**2 = x**2 - 2*x*q_c*c1 + q_c**2*c1**2
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
    return theta1_1
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

t = Symbol('t')
theta2 = dynamicsymbols('theta2')
theta2_d = diff(theta2, t)
theta2_dd = diff(theta2_d, t)
theta3 = dynamicsymbols('theta3')
theta3_d = diff(theta3, t)
theta3_dd = diff(theta3_d, t)
theta4 = dynamicsymbols('theta4')
theta4_d = diff(theta4, t)
theta4_dd = diff(theta4_d, t)

theta1 = findTheta1(theta2, theta3, STEP_POSITION - w_r)
theta1_d = diff(theta1, t)
theta1_dd = diff(theta1_d, t)

theta = Matrix([theta2, theta3, theta4])
theta_d = Matrix([theta2_d, theta3_d, theta4_d])
theta_dd = Matrix([theta2_dd, theta3_dd, theta4_dd])

theta2_0 = Symbol('theta2_0')
theta2_d_0 = Symbol('theta2_d_0')
theta2_dd_0 = Symbol('theta2_dd_0')
theta3_0 = Symbol('theta3_0')
theta3_d_0 = Symbol('theta3_d_0')
theta3_dd_0 = Symbol('theta3_dd_0')
theta4_0 = Symbol('theta4_0')
theta4_d_0 = Symbol('theta4_d_0')
theta4_dd_0 = Symbol('theta4_dd_0')

theta_0 = [theta2_0, theta3_0, theta4_0]
theta_d_0 = [theta2_d_0, theta3_d_0, theta4_d_0]

# ORDER IS IMPORTANT!
# We want to substitute theta2_d before theta2 because Derivative(theta2, t) != Derivative(theta2_0, t)
# Substitution will fail if we substituted theta2 first
# Should no longer be necessary: https://github.com/sympy/sympy/issues/17656
substitutions = [
    (theta2_dd, theta2_dd_0),
    (theta3_dd, theta3_dd_0),
    (theta4_dd, theta4_dd_0),
    (theta2_d, theta2_d_0),
    (theta3_d, theta3_d_0),
    (theta4_d, theta4_d_0),
    (theta2, theta2_0),
    (theta3, theta3_0),
    (theta4, theta4_0)]

theta_dd_0 = [theta2_dd_0, theta3_dd_0, theta4_dd_0]

tau2 = Symbol('tau2')
tau3 = Symbol('tau3')
tau4 = Symbol('tau4')
tau = [tau2, tau3, tau4]

lambda_parameters = [theta_0, theta_d_0, tau]

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
Pb = P2 + (l_2/2.0*c12 + l_b*cos(theta1 + theta2 + pi/2.0))*i - (l_2/2.0*s12 + l_b*sin(theta1 + theta2 + pi/2.0))*k

I_w = 1.0/2.0*m_w*(w_r**2)
# I_1 = m_2*P2.dot(P2) + m_3*P3.dot(P3) + m_b*Pb.dot(Pb) + m_4*P4.dot(P4)
I_2 = m_3*(P3-P2).dot(P3-P2) + m_b*(Pb-P2).dot(Pb-P2) + m_4*(P4-P2).dot(P4-P2)
I_3 = m_4*(P4-P3).dot(P4-P3)
I_4 = I_w
I = Matrix([I_2, I_3, I_4])

'''
Generate calc_end_effector_force()
'''
tic = time.time()
J = Matrix([P4.dot(i), P4.dot(k)]).jacobian(theta)
# https://studywolf.wordpress.com/2013/09/02/robot-control-jacobians-velocity-and-force/
# x_d = J*q_d
# x_dd = J*q_dd + J_d*q_d
# F = m*(J*q_dd + J_d*q_d)
# F = m*(J*(tau/I) + J_d*q_d)
q_dd = Matrix([tau[i] / I[i] for i in range(len(tau))])
end_effector_force = m_4*((J*q_dd) + J.diff(t)*theta_d)
ee_force_np = lambdify(
        lambda_parameters,
        end_effector_force,
        modules=numpy_trig)
ee_force_sym = lambdify(
        lambda_parameters,
        end_effector_force,
        modules=symbolic_trig)

def calc_end_effector_force(
        theta2, theta3, theta4,
        theta2_d, theta3_d, theta4_d,
        tau2, tau3, tau4, is_symbolic = False):
    if is_symbolic:
        ee_force = ee_force_sym
    else:
        ee_force = ee_force_np
    return ee_force(
            (theta2, theta3, theta4),
            (theta2_d, theta3_d, theta4_d),
            (tau2, tau3, tau4))

logDuration("Generate end_effector_from_torques()", tic)

'''
Generate EOM
'''
tic = time.time()

v2 = P2.diff(t, N)
v3 = P3.diff(t, N)
v4 = P4.diff(t, N)
vb = Pb.diff(t, N)

KE = 0.5*m_2*v2.dot(v2) + 0.5*m_3*v3.dot(v3) + 0.5*m_4*v4.dot(v4) + 0.5*m_b*vb.dot(vb)
PE = m_2*g*P2.dot(k) + m_3*g*P3.dot(k) + m_4*g*P4.dot(k) + m_b*g*Pb.dot(k)
L = Matrix([KE - PE])

lhs = (L.jacobian(theta_d).diff(t) - L.jacobian(theta)).T

# External forces
force = tau4*w_r*k # Vertical force input force (from wheel)

rhs1 = tau2 + force.dot(P4.diff(theta2, N))
rhs2 = tau3 + force.dot(P4.diff(theta3, N))
rhs3 = tau4 + force.dot(P4.diff(theta4, N))
rhs = Matrix([rhs1, rhs2, rhs3])

eom = Eq(lhs, rhs).subs(substitutions)
logDuration("Formulate Lagrange", tic)

# tic = time.time()
# eom = trigsimp(eom, method="fu") # Simplify by minimizing trigonometric functions
# logDuration("Simplify", tic)

tic = time.time()
theta_dd_eom = solve(eom, theta_dd_0, dict=True, simplify=False, rational=False)
logDuration("Solve EOM", tic)

tic = time.time()
theta2_dd_eom = theta_dd_eom[0][theta2_dd_0]
theta3_dd_eom = theta_dd_eom[0][theta3_dd_0]
theta4_dd_eom = theta_dd_eom[0][theta4_dd_0]

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

theta4_dd_np = lambdify(
        lambda_parameters,
        theta4_dd_eom,
        modules=numpy_trig)
theta4_dd_sym = lambdify(
        lambda_parameters,
        theta4_dd_eom,
        modules=symbolic_trig)

logDuration("Lambidfy EOM", tic)

def calc_theta2_dd(
        theta2, theta3, theta4,
        theta2_d, theta3_d, theta4_d,
        tau2, tau3, tau4, is_symbolic = False):
    if is_symbolic:
        theta2_dd_lambd = theta2_dd_sym
    else:
        theta2_dd_lambd = theta2_dd_np
    return theta2_dd_lambd(
            (theta2, theta3, theta4),
            (theta2_d, theta3_d, theta4_d),
            (tau2, tau3, tau4))

def calc_theta3_dd(
        theta2, theta3, theta4,
        theta2_d, theta3_d, theta4_d,
        tau2, tau3, tau4, is_symbolic = False):
    if is_symbolic:
        theta3_dd_lambd = theta3_dd_sym
    else:
        theta3_dd_lambd = theta3_dd_np
    return theta3_dd_lambd(
            (theta2, theta3, theta4),
            (theta2_d, theta3_d, theta4_d),
            (tau2, tau3, tau4))

def calc_theta4_dd(
        theta2, theta3, theta4,
        theta2_d, theta3_d, theta4_d,
        tau2, tau3, tau4, is_symbolic = False):
    if is_symbolic:
        theta4_dd_lambd = theta4_dd_sym
    else:
        theta4_dd_lambd = theta4_dd_np
    return theta4_dd_lambd(
            (theta2, theta3, theta4),
            (theta2_d, theta3_d, theta4_d),
            (tau2, tau3, tau4))

if __name__ == "__main__":
    from random import random

    tic = time.time()
    for i in range(100):
        calc_theta2_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random())
        calc_theta3_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random())
        calc_theta4_dd(
                random(), random(), random(),
                random(), random(), random(),
                random(), random(), random())
    print("300 EOM calculation in " + str(time.time() - tic) + "s")

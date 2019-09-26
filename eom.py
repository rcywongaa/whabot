import time

tic = time.time()
from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy.utilities.lambdify import lambdify
from sympy import diff, Symbol, sin, cos, pi, Matrix
from sympy import Eq, solve
from sympy import init_printing, pprint, pretty
# print("Imported sympy in " + str(time.time() - tic) + "s")

import yaml
import pdb

import pydrake.symbolic

custom_trig = [{
        'sin': pydrake.symbolic.sin,
        'cos': pydrake.symbolic.cos}, 'numpy']

with open("constants.yaml", 'r') as stream:
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
    (theta1_d, theta1_d_0),
    (theta2_d, theta2_d_0),
    (theta3_d, theta3_d_0),
    (theta1, theta1_0),
    (theta2, theta2_0),
    (theta3, theta3_0)]

F = Symbol('F')
tau1 = Symbol('tau1')
tau2 = Symbol('tau2')
tau3 = Symbol('tau3')
# tau4 = Symbol('tau4')

m1 = constants['wheel_mass']
m2 = constants['motor1_mass']
m3 = constants['motor2_mass']
m4 = constants['wheel_mass']
mb = constants['battery_mass']
l1 = constants['link1_length']
l2 = constants['link2_length']
l3 = constants['link3_length']
lb = constants['linkb_length']

# m1 = Symbol('m1')
# m2 = Symbol('m2')
# m3 = Symbol('m3')
# m4 = Symbol('m4')
# mb = Symbol('mb')
# l1 = Symbol('l1')
# l2 = Symbol('l2')
# l3 = Symbol('l3')
# lb = Symbol('lb')

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

P2 = l1*c1*i + l1*s1*j
P3 = P2 + l2*c12*i + l2*s12*j
P4 = P3 + l3*c123*i + l2*s123*j
Pb = P2 + l2/2.0*c12*i + l2/2.0*s12*j + lb*cos(theta1 + theta2 + pi/2.0)*i + lb*sin(theta1 + theta2 + pi/2.0)*j

v2 = P2.diff(t, N)
v3 = P3.diff(t, N)
v4 = P4.diff(t, N)
vb = Pb.diff(t, N)

KE = 0.5*m1*v2.dot(v2) + 0.5*m2*v2.dot(v2) + 0.5*m3*v3.dot(v3) + 0.5*m4*v4.dot(v4) + 0.5*mb*vb.dot(vb)
PE = m2*g*P2.dot(j) + m3*g*P3.dot(j) + m4*g*P4.dot(j) + mb*g*Pb.dot(j)
L = KE - PE

lhs = (Matrix([L]).jacobian(theta_d).diff(t) - Matrix([L]).jacobian(theta)).T

rhs1 = F*(l1*s1 + l2*s12 + l3*s123)
rhs2 = tau2 + F*(l2*s12 + l3*s123)
rhs3 = tau3 + F*(l3*s123)
rhs = Matrix([rhs1, rhs2, rhs3])

eom = Eq(lhs, rhs)
# print("Lagrange formulated in " + str(time.time() - tic) + "s")

tic = time.time()
theta_dd_eom = solve(eom, theta_dd, dict=True, simplify=False)
# print("EOM solved in " + str(time.time() - tic) + "s")

tic = time.time()
theta1_dd_eom = theta_dd_eom[0][theta1_dd].subs(substitutions)
theta2_dd_eom = theta_dd_eom[0][theta2_dd].subs(substitutions)
theta3_dd_eom = theta_dd_eom[0][theta3_dd].subs(substitutions)
theta1_dd_lambd = lambdify(
        [theta_0, theta_d_0, (tau1, tau2, tau3, F)],
        theta1_dd_eom,
        modules=custom_trig)
theta2_dd_lambd = lambdify(
        [theta_0, theta_d_0, (tau1, tau2, tau3, F)],
        theta2_dd_eom,
        modules=custom_trig)
theta3_dd_lambd = lambdify(
        [theta_0, theta_d_0, (tau1, tau2, tau3, F)],
        theta3_dd_eom,
        modules=custom_trig)
# print("Lambdified in " + str(time.time() - tic) + "s")

def calc_theta1_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F):
    return theta1_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

def calc_theta2_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F):
    return theta2_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

def calc_theta3_dd(
        theta1, theta2, theta3,
        theta1_d, theta2_d, theta3_d,
        tau1, tau2, tau3, F):
    return theta3_dd_lambd(
            (theta1, theta2, theta3),
            (theta1_d, theta2_d, theta3_d),
            (tau1, tau2, tau3, F))

if __name__ == "__main__":
    tic = time.time()
    from random import random
    for i in range(100):
        theta1_dd_lambd(
                (random(), random(), random()),
                (random(), random(), random()),
                (random(), random(), random(), random()))
        theta2_dd_lambd(
                (random(), random(), random()),
                (random(), random(), random()),
                (random(), random(), random(), random()))
        theta3_dd_lambd(
                (random(), random(), random()),
                (random(), random(), random()),
                (random(), random(), random(), random()))
    print("300 EOM calculation in " + str(time.time() - tic) + "s")
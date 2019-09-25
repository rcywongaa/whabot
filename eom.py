from sympy.physics.vector import dynamicsymbols, ReferenceFrame, dot
from sympy import diff, Symbol, sin, cos, pi, Matrix
from sympy import init_printing, pprint, pretty
import pdb

g = 9.81

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

m1 = Symbol('m1')
m2 = Symbol('m2')
m3 = Symbol('m3')
m4 = Symbol('m4')
mb = Symbol('mb')
l1 = Symbol('l1')
l2 = Symbol('l2')
l3 = Symbol('l3')
lb = Symbol('lb')

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

pdb.set_trace()
lhs = L.diff(theta_d).diff(t) - L.diff(theta)

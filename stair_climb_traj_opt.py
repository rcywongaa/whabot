import pydrake
from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
from pydrake.symbolic import Variable
import numpy as np
import math
from enum import Enum
import time
import pdb

from pydrake.geometry import ConnectDrakeVisualizer, SceneGraph, HalfSpace, Box
from pydrake.lcm import DrakeLcm
from pydrake.multibody.tree import UniformGravityFieldElement
from pydrake.multibody.plant import MultibodyPlant, AddMultibodyPlantSceneGraph, CoulombFriction
from pydrake.multibody.parsing import Parser
from pydrake.systems.framework import DiagramBuilder
from pydrake.systems.analysis import Simulator
from pydrake.math import RigidTransform

import eom

STATE_SIZE = 8
TORQUE_SIZE = 3
STEP_DEPTH = 0.4
STEP_WIDTH = 0.5
STEP_HEIGHT = 0.3
HUB_MOTOR_MAX_TORQUE = 20
MIN_NORMAL_REACTION = 4
DISCRIMINANT_EPSILON = 1e-10
DYNAMICS_EPSILON = 1e-3
g = 9.81
l_1 = 0.8
l_2 = 0.32
l_3 = l_1
l_b = 0.6*l_1 # length of battery link
m_2 = 0.2 # mass of l_1 - l_2 link motor
m_m = 0.4 # mass of battery link motor
m_b = 4.0 # mass of battery
m_3 = m_2 # mass of l_2 - l_3 motor
m_w = 4.0 # mass of wheel
w_r = 0.1 # wheel radius
I_w = 1.0/2.0*m_w*(w_r**2)
NUM_TIME_STEPS = 100
TIME_INTERVAL = 1.0/NUM_TIME_STEPS
# F_f <= COEFF_FRICTION * F_n
COEFF_FRICTION = 0.6

def calc_theta1_dd(state, u):
    return eom.calc_theta1_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def calc_theta2_dd(state, u):
    return eom.calc_theta2_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def calc_theta3_dd(state, u):
    return eom.calc_theta3_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def findTheta1(theta2, theta3, theta4, is_symbolic = True):
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
    q_c = l_2*s2 + l_3*s23
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
    y = l_1*s1 + l_2*s12 + l_3*s123
    return (x, y)

def derivs(state, tau234, is_symbolic = True):
    if is_symbolic:
        sin = pydrake.symbolic.sin
        cos = pydrake.symbolic.cos
    else:
        sin = np.sin
        cos = np.cos

    tau2 = tau234[0]
    tau3 = tau234[1]
    tau4 = tau234[2]
    theta1 = state[0]
    theta2 = state[1]
    theta3 = state[2]
    theta4 = state[3]
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

    tau1 = 0.0 # Unactuated
    tau = np.array([tau1, tau2, tau3, tau4])
    state_d = np.zeros_like(state)
    state_d[0:4] = state[4:8]
    state_d[4] = calc_theta1_dd(state, tau)
    state_d[5] = calc_theta2_dd(state, tau)
    state_d[6] = calc_theta3_dd(state, tau)
    state_d[7] = tau234[2] / I_w
    return state_d

if __name__ == "__main__":
    mp = MathematicalProgram()
    state_over_time = np.zeros(shape=(NUM_TIME_STEPS, STATE_SIZE), dtype=pydrake.symbolic.Variable)
    tau234_over_time = np.zeros(shape=(NUM_TIME_STEPS, TORQUE_SIZE), dtype=pydrake.symbolic.Variable)

    initial_state = mp.NewContinuousVariables(8, "state_0")

    # Constrain initial velocity to be 0
    for j in range(4, 8):
        mp.AddConstraint(initial_state[j] <= 0.0)
        mp.AddConstraint(initial_state[j] >= 0.0)

    # Constrain initial theta to be between 0 ~ 180
    for j in range(0, 4):
        mp.AddConstraint(initial_state[j] >= 0.0)
        mp.AddConstraint(initial_state[j] <= np.pi)

    # for j in range(0, 3):
        # mp.AddConstraint(initial_state[j] <= np.pi/2.0)
        # mp.AddConstraint(initial_state[j] >= np.pi/2.0)

    # Constrain initial theta4 (i.e. front wheel) to be 0.0
    # mp.AddConstraint(initial_state[3] <= 0.0)
    # mp.AddConstraint(initial_state[3] >= 0.0)

    # Constrain initial front wheel position y position to be 0.0
    initial_wheel_position = findFrontWheelPosition(initial_state[0], initial_state[1], initial_state[2])
    mp.AddConstraint(initial_wheel_position[1] <= 0.0)
    mp.AddConstraint(initial_wheel_position[1] >= 0.0)

    state_over_time[0] = initial_state

    for i in range(NUM_TIME_STEPS-1):
        print("Adding constraints for t = " + str(i))

        tau234 = mp.NewContinuousVariables(3, "tau234_%d" % i)
        next_state = mp.NewContinuousVariables(8, "state_%d" % (i+1))

        theta1 = next_state[0]
        theta2 = next_state[1]
        theta3 = next_state[2]
        theta4 = next_state[3]
        theta1_d = next_state[4]
        theta2_d = next_state[5]
        theta3_d = next_state[6]
        theta4_d = next_state[7]
        tau1 = 0.0 # Unactuated
        tau2 = tau234[0]
        tau3 = tau234[1]
        tau4 = tau234[2]

        # Add end force constraint
        # eom.calc_end_force_from_torques(
                # theta1, theta2, theta3,
                # theta1_d, theta2_d, theta3_d,
                # tau1, tau2, tau3)
        # mp.AddConstraint(COEFF_FRICTION*J.dot(tau123)[0] <= -w_r*tau234[2]) # FIXME

        # Constrain no x motion of front wheel
        # mp.AddConstraint(findFrontWheelPosition(next_state[0], next_state[1], next_state[2])[0] <= initial_wheel_position[0])
        # mp.AddConstraint(findFrontWheelPosition(next_state[0], next_state[1], next_state[2])[0] >= initial_wheel_position[0])

        # for j in range(3): # Constrain theta1, theta2, theta3
            # mp.AddConstraint(next_state[j] >= 0.0)
            # mp.AddConstraint(next_state[j] <= np.pi)
        # for j in range(8):
            # mp.AddConstraint(next_state[j] <= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])
            # mp.AddConstraint(next_state[j] >= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])

        def next_state_constraint(stacked):
            next_state = stacked[0:STATE_SIZE]
            current_state = stacked[STATE_SIZE:STATE_SIZE*2]
            current_tau234 = stacked[STATE_SIZE*2:STATE_SIZE*2+TORQUE_SIZE]
            desired_state = current_state + TIME_INTERVAL*derivs(current_state, current_tau234, is_symbolic = False)
            diff = desired_state - next_state
            ret = np.zeros_like(stacked)
            ret[0:STATE_SIZE] = diff
            return ret

        stacked = np.concatenate([next_state, state_over_time[i], tau234])
        bounds = np.ones(stacked.shape)*DYNAMICS_EPSILON
        mp.AddConstraint(next_state_constraint, -bounds, bounds, stacked)

        tau234_over_time[i] = tau234
        state_over_time[i+1] = next_state

    mp.AddCost(0.01 * tau234_over_time[:,0].dot(tau234_over_time[:,0]))
    mp.AddCost(0.01 * tau234_over_time[:,1].dot(tau234_over_time[:,1]))
    mp.AddCost(0.01 * tau234_over_time[:,2].dot(tau234_over_time[:,2]))
    # This cost is incorrect, different case for if front wheel is left / right of back wheel
    # mp.AddCost(-(state_over_time[:,0].dot(state_over_time[:,0])))

    final_state = state_over_time[-1]
    # Constrain final velocity to be 0
    for j in range(4, 8):
        mp.AddConstraint(final_state[j] <= 0.0)
        mp.AddConstraint(final_state[j] >= 0.0)

    # Constrain final front wheel position
    final_front_wheel_pos = findFrontWheelPosition(final_state[0], final_state[1], final_state[2])
    mp.AddConstraint(final_front_wheel_pos[1] <= STEP_HEIGHT)
    mp.AddConstraint(final_front_wheel_pos[1] >= STEP_HEIGHT)

    print("Begin solving...")
    t = time.time()
    result = Solve(mp)
    print("Done solving in " + str(time.time() - t) + "s")
    is_success = result.is_success()
    print("is_success = " + str(is_success))
    torque_over_time_star = result.GetSolution(tau234_over_time)
    state_over_time_star = result.GetSolution(state_over_time)
    pdb.set_trace()





    file_name = "res/stair_climb.sdf"
    builder = DiagramBuilder()
    stair_climb, scene_graph = AddMultibodyPlantSceneGraph(builder)
    # stair_climb.RegisterAsSourceForSceneGraph(scene_graph)
    Parser(plant=stair_climb).AddModelFromFile(file_name)
    # stair_climb.AddForceElement(UniformGravityFieldElement())

    initial_state = state_over_time[0]
    front_wheel_x, front_wheel_y = findFrontWheelPosition(initial_state[0], initial_state[1], initial_state[2])
    front_wheel_x = front_wheel_x.Evaluate()
    front_wheel_y = front_wheel_y.Evaluate()
    step = Box(STEP_DEPTH, STEP_WIDTH, STEP_HEIGHT)
    step_pos = RigidTransform([front_wheel_x + w_r + STEP_DEPTH/2.0, 0.0, STEP_HEIGHT/2.0])

    stair_climb.RegisterCollisionGeometry(
            stair_climb.world_body(),
            RigidTransform([0.0, 0.0, 0.0]),
            HalfSpace(),
            "GroundCollision",
            CoulombFriction(COEFF_FRICTION, COEFF_FRICTION))

    stair_climb.RegisterVisualGeometry(
            stair_climb.world_body(),
            RigidTransform([0.0, 0.0, 0.0]),
            HalfSpace(),
            "GroundVisual",
            np.array([0.5, 0.5, 0.5, 0.5])) # Color

    stair_climb.RegisterCollisionGeometry(
            stair_climb.world_body(),
            step_pos,
            step,
            "StepCollision",
            CoulombFriction(COEFF_FRICTION, COEFF_FRICTION))

    stair_climb.RegisterVisualGeometry(
            stair_climb.world_body(),
            step_pos,
            step,
            "StepVisual",
            np.array([1.0, 1.0, 0.0, 1.0])) # Color

    stair_climb.Finalize()

    ConnectDrakeVisualizer(builder=builder, scene_graph=scene_graph)
    diagram = builder.Build()

    diagram_context = diagram.CreateDefaultContext()
    stair_climb_context = diagram.GetMutableSubsystemContext(stair_climb, diagram_context)

    stair_climb_context.FixInputPort(stair_climb.get_actuation_input_port().get_index(), [0, 0, 0, 0, 0])

    theta1 = stair_climb.GetJointByName("theta1")
    theta2 = stair_climb.GetJointByName("theta2")
    theta3 = stair_climb.GetJointByName("theta3")
    theta4 = stair_climb.GetJointByName("theta4")
    phi = stair_climb.GetJointByName("phi")
    theta1.set_angle(context=stair_climb_context, angle=initial_state[0])
    theta2.set_angle(context=stair_climb_context, angle=initial_state[1])
    theta3.set_angle(context=stair_climb_context, angle=initial_state[2])
    theta4.set_angle(context=stair_climb_context, angle=initial_state[3])
    phi.set_angle(context=stair_climb_context, angle=initial_state[4])

    simulator = Simulator(diagram, diagram_context)
    simulator.set_publish_every_time_step(False)
    simulator.set_target_realtime_rate(0.1)
    simulator.Initialize()
    simulator.AdvanceTo(1.0)

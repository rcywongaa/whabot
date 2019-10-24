import pydrake
from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve, GetInfeasibleConstraints
from pydrake.symbolic import Variable
import numpy as np
import math
from enum import Enum
import time
import pdb
from constants import *

import eom
import visualize

use_symbolic_derivs = False
STATE_SIZE = 8
TORQUE_SIZE = 3

JOINT_MOTOR_MAX_TORQUE = 2.0
HUB_MOTOR_MAX_TORQUE = 20
MIN_NORMAL_REACTION = 4
DYNAMICS_EPSILON = 1e-5

TIME_ALLOWED = 1.0
NUM_TIME_STEPS = 1000
TIME_INTERVAL = TIME_ALLOWED / NUM_TIME_STEPS

def calc_theta2_dd(state, u, is_symbolic):
    return eom.calc_theta2_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2],
            is_symbolic)

def calc_theta3_dd(state, u, is_symbolic):
    return eom.calc_theta3_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2],
            is_symbolic)

def calc_theta4_dd(state, u, is_symbolic):
    return eom.calc_theta1_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2],
            is_symbolic)

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

    tau = np.array([tau2, tau3, tau4])
    state_d = np.zeros_like(state)
    state_d[0:4] = state[4:8]
    state_d[4] = calc_theta1_dd(state, tau, is_symbolic)
    state_d[5] = calc_theta2_dd(state, tau, is_symbolic)
    state_d[6] = calc_theta3_dd(state, tau, is_symbolic)
    state_d[7] = tau234[2] / I_w
    return state_d

def constrain_theta123(mp, theta1, theta2, theta3):
    # note theta1 should be in opposite direction
    mp.AddConstraint(theta1 <= 0.0).evaluator().set_description("Constrain theta1 <= 0.0")
    mp.AddConstraint(theta1 >= -np.pi).evaluator().set_description("Constrain theta1 >= -np.pi")
    mp.AddConstraint(theta2 >= 0.0).evaluator().set_description("Constrain theta2 >= 0.0")
    mp.AddConstraint(theta2 <= np.pi).evaluator().set_description("Constrain theta2 <= np.pi")
    mp.AddConstraint(theta3 >= 0.0).evaluator().set_description("Constrain theta3 >= 0.0")
    mp.AddConstraint(theta3 <= np.pi).evaluator().set_description("Constrain theta3 <= np.pi")

if __name__ == "__main__":
    mp = MathematicalProgram()
    state_over_time = np.zeros(shape=(NUM_TIME_STEPS, STATE_SIZE), dtype=pydrake.symbolic.Variable)
    tau234_over_time = np.zeros(shape=(NUM_TIME_STEPS, TORQUE_SIZE), dtype=pydrake.symbolic.Variable)

    initial_state = mp.NewContinuousVariables(8, "state_0")
    initial_theta1 = initial_state[0]
    initial_theta2 = initial_state[1]
    initial_theta3 = initial_state[2]
    initial_theta4 = initial_state[3]

    # Constrain initial velocity to be 0
    # for j in range(4, 8):
        # mp.AddConstraint(initial_state[j] <= 0.0).evaluator().set_description("Constrain initial_state[%d] <= 0.0" % j)
        # mp.AddConstraint(initial_state[j] >= 0.0).evaluator().set_description("Constrain initial_state[%d] >= 0.0" % j)

    # constrain_theta123(mp, initial_theta1, initial_theta2, initial_theta3)
    mp.AddConstraint(initial_theta1 <= -np.pi/3.0)
    mp.AddConstraint(initial_theta1 >= -np.pi/3.0)
    mp.AddConstraint(initial_theta2 <= np.pi/3.0)
    mp.AddConstraint(initial_theta2 >= np.pi/3.0)
    # mp.AddConstraint(initial_theta4 >= 0.0)
    # mp.AddConstraint(initial_theta4 <= 0.0)

    # Constrain initial front wheel position y position to be 0.0
    initial_wheel_position = eom.findFrontWheelPosition(initial_state[0], initial_state[1], initial_state[2])
    mp.AddConstraint(initial_wheel_position[1] <= 0.0).evaluator().set_description("Constrain initial_wheel_position[1] <= 0.0")
    mp.AddConstraint(initial_wheel_position[1] >= 0.0).evaluator().set_description("Constrain initial_wheel_position[1] >= 0.0")

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
        tau2 = tau234[0] # joint motor
        tau3 = tau234[1] # joint motor
        tau4 = tau234[2] # hub motor

        # Add torque constraints
        # mp.AddConstraint(tau2**2 <= JOINT_MOTOR_MAX_TORQUE**2)
        # mp.AddConstraint(tau3**2 <= JOINT_MOTOR_MAX_TORQUE**2)
        # mp.AddConstraint(tau4**2 <= HUB_MOTOR_MAX_TORQUE**2)

        # constrain_theta123(mp, theta1, theta2, theta3)

        if use_symbolic_derivs:
            for j in range(8):
                mp.AddConstraint(next_state[j] <= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])
                mp.AddConstraint(next_state[j] >= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])
        else:
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
            mp.AddConstraint(next_state_constraint, -bounds, bounds, stacked).evaluator().set_description("Constrain next_state_constraint %d" % i)

        tau234_over_time[i] = tau234
        state_over_time[i+1] = next_state

    mp.AddCost(0.1 * tau234_over_time[:,0].dot(tau234_over_time[:,0]))
    mp.AddCost(0.1 * tau234_over_time[:,1].dot(tau234_over_time[:,1]))
    mp.AddCost(0.1 * tau234_over_time[:,2].dot(tau234_over_time[:,2]))
    # This cost is incorrect, different case for if front wheel is left / right of back wheel
    # mp.AddCost(-(state_over_time[:,0].dot(state_over_time[:,0])))

    final_state = state_over_time[-1]
    # Constrain final velocity to be 0
    # for j in range(4, 8):
        # mp.AddConstraint(final_state[j] <= 0.0)
        # mp.AddConstraint(final_state[j] >= 0.0)

    # Constrain final front wheel position
    final_front_wheel_pos = eom.findFrontWheelPosition(final_state[0], final_state[1], final_state[2])
    # mp.AddConstraint(final_front_wheel_pos[1] <= STEP_HEIGHT - w_r)
    mp.AddConstraint(final_front_wheel_pos[1] >= STEP_HEIGHT - w_r).evaluator().set_description("Constrain final_front_wheel_pos[1] >= 0.01")

    print("Begin solving...")
    t = time.time()
    result = Solve(mp)
    solve_traj_opt_duration = time.time() - t
    print("Done solving in " + str(solve_traj_opt_duration) + "s (" + str(solve_traj_opt_duration/60.0) + "m)")
    is_success = result.is_success()
    print("is_success = " + str(is_success))
    torque_over_time_star = result.GetSolution(tau234_over_time)
    state_over_time_star = result.GetSolution(state_over_time)

    time_step = 0
    last_input = 'c'
    max_time_step = len(state_over_time_star)
    while time_step < max_time_step:
        state = state_over_time_star[time_step]
        torque = np.array([i.Evaluate() for i in torque_over_time_star[time_step]])
        state_d = derivs(state, torque, is_symbolic = False)
        proposed_next_state = state + TIME_INTERVAL*state_d
        print("time_step = " + str(time_step))
        print("theta = " + str(state[0:int(STATE_SIZE/2)]))
        print("theta_d = " + str(state[int(STATE_SIZE/2):STATE_SIZE]))
        print("torques = " + str(torque))
        print("theta_dd = " + str(state_d[int(STATE_SIZE/2):STATE_SIZE]))
        print("proposed_next_state = " + str(proposed_next_state))
        visualize.visualize(state[0], state[1], state[2])
        print("c = continue, r = reverse, q = quit")
        user_input = input("Input: ")

        if user_input == '':
            user_input = last_input
        if user_input == 'c':
            time_step += 1
        elif user_input == 'r':
            time_step = max(0, time_step - 1)
        elif user_input == 'q':
            break
        else:
            try:
                increment = int(user_input)
            except ValueError:
                print("Not a number")
                continue
            time_step = max(min(time_step + increment, max_time_step-1), 0)

        last_input = user_input
    if not is_success:
        # https://drake.mit.edu/pydrake/pydrake.solvers.mathematicalprogram.html?highlight=mathematicalprogram#pydrake.solvers.mathematicalprogram.GetInfeasibleConstraints
        infeasible_constraints = GetInfeasibleConstraints(mp, result)
        # https://drake.mit.edu/pydrake/pydrake.solvers.mathematicalprogram.html?highlight=mathematicalprogram#pydrake.solvers.mathematicalprogram.MathematicalProgram.AddVisualizationCallback
        # https://drake.mit.edu/pydrake/pydrake.solvers.mathematicalprogram.html?highlight=mathematicalprogram#pydrake.solvers.mathematicalprogram.MathematicalProgramResult.GetSuboptimalSolution
    pdb.set_trace()

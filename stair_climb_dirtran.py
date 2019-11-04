import pydrake
from pydrake.multibody.plant import MultibodyPlant, AddMultibodyPlantSceneGraph, CoulombFriction
from pydrake.multibody.tree import UniformGravityFieldElement
from pydrake.multibody.tree import ForceElement
from pydrake.multibody.parsing import Parser
from pydrake.systems.framework import DiagramBuilder
from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve, GetInfeasibleConstraints
from pydrake.symbolic import Variable
import numpy as np
import math
from enum import Enum
import time
import pdb
from constants import *

import constraint_equation
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

def constrain_theta123(mp, theta1, theta2, theta3):
    # note theta1 should be in opposite direction
    mp.AddConstraint(theta1 <= 0.0).evaluator().set_description("Constrain theta1 <= 0.0")
    mp.AddConstraint(theta1 >= -np.pi).evaluator().set_description("Constrain theta1 >= -np.pi")
    mp.AddConstraint(theta2 >= 0.0).evaluator().set_description("Constrain theta2 >= 0.0")
    mp.AddConstraint(theta2 <= np.pi).evaluator().set_description("Constrain theta2 <= np.pi")
    mp.AddConstraint(theta3 >= 0.0).evaluator().set_description("Constrain theta3 >= 0.0")
    mp.AddConstraint(theta3 <= np.pi).evaluator().set_description("Constrain theta3 <= np.pi")

# class FrontWheelForce(ForceElement):
    # def __init__(self, plant):
        # front_wheel = plant.GetBodyByName("front_wheel")
        # self.input_port_index = 0
        # front_wheel_node_index = front_wheel.index()
        # pdb.set_trace()
        # ForceElement.__init__(self, front_wheel.model_instance())

    # def calcForce(self, context):
        # current_tau4 = self.EvalVectorInput(context, self.input_port_index)[4]
        # return w_r*tau4

    # def DoCalcAndAddForceContribution(self, context, pc, vc, forces):
        # f = self.calcForce(context)
        # mutable_forces = forces.mutable_body_forces()
        # spatial_force = SpatialForce([0, 0, 0], [0, 0, f])
        # mutable_forces[self.front_wheel_node_index] += spatial_force

    # def CalcPotentialEnergy(self, context, pc):
        # return 0.0

    # def CalcConservativePower(self, context, pc, vc):
        # return 0.0

    # def CalcNonConservativePower(self, context, pc, vc):
        # vel = vc.get_V_WB(self.front_wheel_node_index)
        # z_d = vel.translational()[2]
        # return self.calcForce(context)*z_d

class ConstraintForce(ForceElement):
    def __init__(self, plant):
        self.plant = plant
        self.theta1_joint_index = plant.GetJointByName("theta1").index()
        ForceElement.__init__(self, plant.GetBodyByName("front_wheel").model_instance())

    def calcTorque(self, context):
        pdb.set_trace()
        q = plant.GetPositions(context)
        q_d = plant.GetVelocities(context)
        H = constraint_equation.calc_H(q)
        H_d = constraint_equation.calc_H_d(q, q_d)
        mass_matrix = plant.CalcMassMatrixViaInverseDynamics(context)
        tau_g = plant.CalcGravityGeneralizedForces(context)
        torque = -(H*M.inv()*H.T).pinv()*(H*M.inv()*tau_g + H*q_d)
        return torque

    def DoCalcAndAddForceContribution(self, context, pc, vc, forces):
        tau = self.calcTorque(context)
        mutable_forces = forces.mutable_body_forces()
        spatial_force = SpatialForce([0, 0, tau], [0, 0, 0])
        mutable_forces[self.theta1_joint_index] += spatial_force

    def CalcPotentialEnergy(self, context, pc):
        return 0.0

    def CalcConservativePower(self, context, pc, vc):
        return 0.0

    def CalcNonConservativePower(self, context, pc, vc):
        pass

if __name__ == "__main__":
    file_name = "res/stair_climb.sdf"
    builder = DiagramBuilder()
    stair_climb, scene_graph = AddMultibodyPlantSceneGraph(builder)
    Parser(plant=stair_climb).AddModelFromFile(file_name)
    # stair_climb.AddForceElement(FrontWheelForce(stair_climb))
    stair_climb.AddForceElement(ConstraintForce(stair_climb))
    stair_climb.Finalize()

    context = stair_climb.CreateDefaultContext()

    dirtran = DirectTranscription(stair_climb, context, num_time_samples=100)
    dirtran.AddConstraint(dirtran.initial_state() == [-np.pi/3, np.pi/3, np.pi/3, 0.0, 0.0])
    u = dirtran.input()
    dirtran.AddRunningCost(u.dot(u))
    # dirtran.AddConstraintToAllKnotPoints(u[0] == 0)
    final_state = dirtran.final_state()
    final_front_wheel_position = eom.FindFrontWheelPosition(final_state[0], final_state[1], final_state[2])
    dirtran.AddConstraint(final_front_wheel_position[1] == STEP_HEIGHT - w_r)

    result = Solve(dirtran)

    times = dirtran.GetSampleTimes(result)
    inputs = dirtran.GetInputSamples(result)
    states = dirtran.GetStateSamples(result)
    input_traj = dirtran.ReconstructInputTrajectory(result)
    state_traj = dirtran.ReconstructStateTrajectory(result)
    pdb.set_trace()

    # torque_over_time_star = result.GetSolution(tau234_over_time)
    # state_over_time_star = result.GetSolution(state_over_time)

    # time_step = 0
    # last_input = 'c'
    # max_time_step = len(state_over_time_star)
    # while time_step < max_time_step:
        # state = state_over_time_star[time_step]
        # torque = np.array([i.Evaluate() for i in torque_over_time_star[time_step]])
        # state_d = derivs(state, torque, is_symbolic = False)
        # proposed_next_state = state + TIME_INTERVAL*state_d
        # print("time_step = " + str(time_step))
        # print("theta = " + str(state[0:int(STATE_SIZE/2)]))
        # print("theta_d = " + str(state[int(STATE_SIZE/2):STATE_SIZE]))
        # print("torques = " + str(torque))
        # print("theta_dd = " + str(state_d[int(STATE_SIZE/2):STATE_SIZE]))
        # print("proposed_next_state = " + str(proposed_next_state))
        # visualize.visualize(state[0], state[1], state[2])
        # print("c = continue, r = reverse, q = quit")
        # user_input = input("Input: ")

        # if user_input == '':
            # user_input = last_input
        # if user_input == 'c':
            # time_step += 1
        # elif user_input == 'r':
            # time_step = max(0, time_step - 1)
        # elif user_input == 'q':
            # break
        # else:
            # try:
                # increment = int(user_input)
            # except ValueError:
                # print("Not a number")
                # continue
            # time_step = max(min(time_step + increment, max_time_step-1), 0)

        # last_input = user_input
    # pdb.set_trace()

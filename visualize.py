import numpy as np
from pydrake.geometry import ConnectDrakeVisualizer, SceneGraph, HalfSpace, Box
from pydrake.lcm import DrakeLcm
from pydrake.multibody.tree import UniformGravityFieldElement
from pydrake.multibody.plant import MultibodyPlant, AddMultibodyPlantSceneGraph, CoulombFriction
from pydrake.multibody.parsing import Parser
from pydrake.systems.framework import DiagramBuilder
from pydrake.systems.analysis import Simulator
from pydrake.math import RigidTransform

from constants import *
import eom

def visualize(theta1, theta2, theta3, theta4=0.0, phi=0.0):
    file_name = "res/stair_climb.sdf"
    builder = DiagramBuilder()
    stair_climb, scene_graph = AddMultibodyPlantSceneGraph(builder)
    # stair_climb.RegisterAsSourceForSceneGraph(scene_graph)
    Parser(plant=stair_climb).AddModelFromFile(file_name)
    # stair_climb.AddForceElement(UniformGravityFieldElement())

    front_wheel_x, front_wheel_y = eom.findFrontWheelPosition(theta1, theta2, theta3)
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

    theta1_joint = stair_climb.GetJointByName("theta1")
    theta2_joint = stair_climb.GetJointByName("theta2")
    theta3_joint = stair_climb.GetJointByName("theta3")
    theta4_joint = stair_climb.GetJointByName("theta4")
    phi_joint = stair_climb.GetJointByName("phi")
    theta1_joint.set_angle(context=stair_climb_context, angle=theta1)
    theta2_joint.set_angle(context=stair_climb_context, angle=theta2)
    theta3_joint.set_angle(context=stair_climb_context, angle=theta3)
    theta4_joint.set_angle(context=stair_climb_context, angle=theta4)
    phi_joint.set_angle(context=stair_climb_context, angle=phi)

    simulator = Simulator(diagram, diagram_context)
    simulator.set_publish_every_time_step(False)
    simulator.set_target_realtime_rate(0.001)
    simulator.Initialize()
    simulator.AdvanceTo(0.001)

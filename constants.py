import yaml

with open("res/constants.yaml", 'r') as stream:
    constants = yaml.safe_load(stream)

g = constants['g']
m_w = constants['wheel_mass']
m_1 = m_w
m_2 = constants['motor1_mass']
m_3 = constants['motor2_mass']
m_4 = m_w
m_b = constants['battery_mass'] # Note we treat battery motor and battery mass as a single mass
l_1 = constants['link1_length']
l_2 = constants['link2_length']
l_3 = constants['link3_length']
l_b = constants['linkb_length']
w_r = constants['wheel_radius']
STEP_DEPTH = constants['step_depth']
STEP_WIDTH = constants['step_width']
STEP_HEIGHT = constants['step_height']
COEFF_FRICTION = constants['coeff_friction']

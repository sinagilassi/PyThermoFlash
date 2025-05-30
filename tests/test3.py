# import libs
import pyThermoFlash as ptf
from rich import print
import os
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# version
print(ptf.__version__)
# check version
print(ptdb.__version__)
# check version
print(ptdblink.__version__)

# SECTION: examples
# !Example 10-1, Page 438, Fundamental of Chemical Engineering Thermodynamics, Kevin D. Dahm

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: thermodb directory
thermodb_dir = os.path.join(os.getcwd(), 'tests', 'thermodb')

# ! benzene (C6H6)
C6H6_file = os.path.join(thermodb_dir, 'benzene-1.pkl')
C6H6 = ptdb.load_thermodb(C6H6_file)
# check
print(C6H6.check())

# general-data
res_ = C6H6.select('general-data')
print(type(res_))
print(res_.table_columns)
print(res_.prop_data)

# ! toluene (C7H8)
C7H8_file = os.path.join(thermodb_dir, 'toluene-1.pkl')
C7H8 = ptdb.load_thermodb(C7H8_file)
# check
print(C7H8.check())

# general-data
res_ = C7H8.select('general')
print(type(res_))
print(res_.table_columns)
print(res_.prop_data)

# =======================================
# ! THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('benzene', C6H6)
thub1.add_thermodb('toluene', C7H8)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'tests', 'thermodb_config_link.yml')

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! THERMOFLASH CALCULATION
# =======================================
# SECTION: vle model
# components
components = ['benzene', 'toluene']

# model source
model_source = {
    'datasource': datasource,
    'equationsource': equationsource
}

# init vle
vle = ptf.vle(components, model_source=model_source)
print(type(vle))

# NOTE: check sources
print(vle.datasource)
print(vle.equationsource)

# SECTION: bubble-point pressure calculation
# alpha (Binary Interaction Parameter)
alpha = [
    [0, 0.3],
    [0.3, 0]
]
# a_ij
a_ij = [
    [0.0, -2.885],
    [2.191, 0.0]
]
# b_ij
b_ij = [
    [0.0, 1124],
    [-863.7, 0.0]
]
# c_ij
c_ij = [
    [0.0, 0.0],
    [0.0, 0.0]
]
# d_ij
d_ij = [
    [0.0, 0.0],
    [0.0, 0.0]
]

# activity inputs
activity_inputs = {
    'alpha': alpha,
    'a_ij': a_ij,
    'b_ij': b_ij,
    'c_ij': c_ij,
    'd_ij': d_ij
}

# calculated activity coefficients (for bubble-pressure calculation)
activity_coefficients = {
    'benzene': 1.0,
    'toluene': 1.1
}

# inputs
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'temperature': [80, 'C'],
}

# SECTION: bubble point pressure calculation
# NOTE: raoult's law
res_bp = vle.bubble_pressure(
    inputs=inputs, equilibrium_model='raoult',)
print(res_bp)

# NOTE: modified raoult's law
res_bp = vle.bubble_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_inputs=activity_inputs)
print(res_bp)

res_bp = vle.bubble_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_coefficients=activity_coefficients)
print(res_bp)

# SECTION: dew point pressure calculation
# NOTE: raoult's law
res_dp = vle.dew_pressure(
    inputs=inputs, equilibrium_model='raoult')
print(res_dp)

# NOTE: modified raoult's law
res_dp = vle.dew_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_inputs=activity_inputs)
print(res_dp)

res_dp = vle.dew_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_coefficients=activity_coefficients)
print(res_dp)

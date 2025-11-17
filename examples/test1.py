# import libs
import os
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB import (
    build_component_model_source,
    build_components_model_source,
    build_model_source
)
from pyThermoLinkDB.models import ComponentModelSource, ModelSource
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoDB import build_component_thermodb_from_reference, ComponentThermoDB
from rich import print
# thermo flash
import pyThermoFlash as ptf

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# SECTION: examples
# !Example 10-1, Page 438, Fundamental of Chemical Engineering Thermodynamics, Kevin D. Dahm

# =======================================
# #️⃣ BUILD COMPONENT THERMODB
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: components


# =======================================
# #️⃣ THERMODB LINK CONFIGURATION
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
# #️⃣ THERMOFLASH CALCULATION
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

# SECTION: bubble-temperature point calculation
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
# activity model
activity_inputs = {
    'alpha': alpha,
    'a_ij': a_ij,
    'b_ij': b_ij,
    'c_ij': c_ij,
    'd_ij': d_ij
}

# inputs
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'pressure': [101.3, 'kPa'],
}

# SECTION: bubble-point temperature calculation
# NOTE: raoult's law
res_bp = vle.bubble_temperature(
    inputs=inputs,
    equilibrium_model='raoult')
print(res_bp)

# NOTE: modified raoult's law
res_bp = vle.bubble_temperature(
    inputs=inputs,
    equilibrium_model='modified-raoult',
    activity_model='NRTL',
    activity_inputs=activity_inputs)
print(res_bp)

# SECTION: dew-point temperature calculation
# NOTE: raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='raoult',
    solver_method='fsolve')
print(res_dp)

# NOTE: modified raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='modified-raoult',
    activity_model='NRTL',
    activity_inputs=activity_inputs,
    solver_method='least-squares')
print(res_dp)

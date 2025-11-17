# import libs
import os
from typing import List
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
from model_source import datasource, equationsource

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# #️⃣ THERMOFLASH CALCULATION
# =======================================
# SECTION: vle model
# components
components = ['benzene-l', 'toluene-l']

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
    'benzene-l': 1.0,
    'toluene-l': 1.1
}

# inputs
inputs = {
    'mole_fraction': {'benzene-l': 0.26, 'toluene-l': 0.74},
    'temperature': [80, 'C'],
}

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

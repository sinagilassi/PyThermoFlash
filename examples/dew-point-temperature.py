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
from pyThermoFlash.core import calc_dew_point_temperature
from model_source import model_source_, datasource, equationsource

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# #️⃣ THERMOFLASH CALCULATION
# =======================================
# SECTION: vle model
# NOTE: components
components = ['benzene-l', 'toluene-l']

# NOTE: model source
model_source = {
    'datasource': datasource,
    'equationsource': equationsource
}

# NOTE: init vle
vle = ptf.vle(
    components=components,
    model_source=model_source
)
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
# NOTE: activity model
activity_inputs = {
    'alpha': alpha,
    'a_ij': a_ij,
    'b_ij': b_ij,
    'c_ij': c_ij,
    'd_ij': d_ij
}

# NOTE: inputs
inputs = {
    'mole_fraction': {'benzene-l': 0.26, 'toluene-l': 0.74},
    'pressure': [101.3, 'kPa'],
}

pressure = Pressure(
    value=101.3,
    unit='kPa'
)
benzene = Component(
    name='benzene',
    formula='C6H6',
    state='l',
    mole_fraction=0.26
)
toluene = Component(
    name='toluene',
    formula='C7H8',
    state='l',
    mole_fraction=0.74
)

# SECTION: dew-point temperature calculation
# NOTE: raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='raoult',
    solver_method='fsolve'
)
print(res_dp)

# ! new method
res_dp = calc_dew_point_temperature(
    components=[benzene, toluene],
    pressure=pressure,
    model_source=model_source_,
)
print(res_dp)

# NOTE: modified raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='modified-raoult',
    activity_model='NRTL',
    activity_inputs=activity_inputs,
    solver_method='least-squares')
print(res_dp)

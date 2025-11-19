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
from pyThermoFlash.core import calc_isothermal_flash, is_flashable
from model_source import model_source_, datasource, equationsource

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# #️⃣ THERMOFLASH CALCULATION
# =======================================
# SECTION: vle model
# components
components = ['water-l', 'ethanol-l']

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

# alpha (Binary Interaction Parameter)
alpha = [
    [0, 0.3],
    [0.3, 0]
]
# a_ij
a_ij = [
    [0.0, 3.458],
    [-0.801, 0.0]
]
# b_ij
b_ij = [
    [0.0, -586.1],
    [246.2, 0.0]
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
    'mole_fraction': {'water-l': 0.50, 'ethanol-l': 0.50},
    'temperature': [30.0, 'C'],
    'pressure': [7.0, 'kPa']
}

temperature = Temperature(value=30.0, unit='C')
pressure = Pressure(value=7.0, unit='kPa')
water = Component(
    name='water',
    formula='H2O',
    state='l',
    mole_fraction=0.50
)
ethanol = Component(
    name='ethanol',
    formula='C2H6O',
    state='l',
    mole_fraction=0.50
)

# SECTION: check flashability
res_check = vle.is_flashable(inputs=inputs)
print(res_check)

res_check = is_flashable(
    components=[water, ethanol],
    temperature=temperature,
    pressure=pressure,
    model_source=model_source_,
)
print(res_check)

# SECTION: flash calculation
# NOTE: flash calculation (least_squares)
res_bp = vle.flash_isothermal(inputs=inputs)
print(res_bp)

# ! new method
res_bp = calc_isothermal_flash(
    components=[water, ethanol],
    temperature=temperature,
    pressure=pressure,
    model_source=model_source_,
)
print(res_bp)

# NOTE: flash calculation (minimize): not recommended
res_bp = vle.flash_isothermal(inputs=inputs, solver_method='minimize')
print(res_bp)

# REVIEW: modified raoult's law (soon)
# res_bp = vle.flash_isothermal(
#     inputs=inputs,
#     equilibrium_model='modified-raoult',
#     activity_model='NRTL',
#     solver_method='least_squares',
#     activity_inputs=activity_inputs)
# print(res_bp)

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
from .reference_content_0 import REFERENCE_CONTENT
from .model_source import datasource, equationsource

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# SECTION: examples
# !Example 10-1, Page 438, Fundamental of Chemical Engineering Thermodynamics, Kevin D. Dahm

# # =======================================
# # #️⃣ CONFIGURATION
# # =======================================
# # NOTE: parent directory
# parent_dir = os.path.dirname(os.path.abspath(__file__))
# print(parent_dir)

# # =======================================
# # #️⃣ BUILD COMPONENT
# # =======================================
# # NOTE: components (normal state)
# # benzene
# benzene = Component(
#     name='benzene',
#     formula='C6H6',
#     state='l'
# )

# # toluene
# toluene = Component(
#     name='toluene',
#     formula='C7H8',
#     state='l'
# )

# # =======================================
# # #️⃣ BUILD COMPONENT THERMODB
# # =======================================
# # NOTE: ignore state properties
# ignore_state_props = ['MW', 'VaPr', 'Cp_IG']

# # SECTION: build component thermodb
# # ! benzene
# benzene_thermodb_component_: ComponentThermoDB | None = build_component_thermodb_from_reference(
#     component_name=benzene.name,
#     component_formula=benzene.formula,
#     component_state=benzene.state,
#     reference_content=REFERENCE_CONTENT,
#     ignore_state_props=ignore_state_props,
# )
# # check
# if benzene_thermodb_component_ is None:
#     raise ValueError("Benzene ThermoDB component build failed.")

# # ! toluene
# toluene_thermodb_component_: ComponentThermoDB | None = build_component_thermodb_from_reference(
#     component_name=toluene.name,
#     component_formula=toluene.formula,
#     component_state=toluene.state,
#     reference_content=REFERENCE_CONTENT,
#     ignore_state_props=ignore_state_props,
# )
# # check
# if toluene_thermodb_component_ is None:
#     raise ValueError("Toluene ThermoDB component build failed.")

# # SECTION: build model source
# # NOTE: with partially matched
# model_source_list: List[ComponentModelSource] = build_components_model_source(
#     components_thermodb=[
#         benzene_thermodb_component_,
#         toluene_thermodb_component_
#     ]
# )

# # NOTE: build model source
# model_source_: ModelSource = build_model_source(
#     source=model_source_list
# )

# # =======================================
# # #️⃣ THERMODB LINK CONFIGURATION
# # =======================================
# # NOTE: build datasource & equationsource
# datasource = model_source_.data_source
# equationsource = model_source_.equation_source

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

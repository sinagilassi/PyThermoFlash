# import libs
import logging
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
from reference_content import REFERENCE_CONTENT

# NOTE: set logging level
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# version
print(ptf.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# #️⃣ CONFIGURATION
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# =======================================
# #️⃣ BUILD COMPONENT
# =======================================
# NOTE: components (normal state)
# benzene
benzene = Component(
    name='benzene',
    formula='C6H6',
    state='l'
)

# toluene
toluene = Component(
    name='toluene',
    formula='C7H8',
    state='l'
)

# ethanol
ethanol = Component(
    name='ethanol',
    formula='C2H6O',
    state='l'
)

# water
water = Component(
    name='water',
    formula='H2O',
    state='l'
)

# components list
components_list: List[Component] = [
    ethanol,
    water
]

# selected component
selected_component: Component = ethanol

# =======================================
# #️⃣ BUILD COMPONENT THERMODB
# =======================================
# NOTE: ignore state properties
ignore_state_props = ['MW', 'VaPr', 'Cp_IG']

# SECTION: build component thermodb
thermodb_component_ = build_component_thermodb_from_reference(
    component_name=selected_component.name,
    component_formula=selected_component.formula,
    component_state=selected_component.state,
    reference_content=REFERENCE_CONTENT,
    ignore_state_props=ignore_state_props,
)
# >> check
if thermodb_component_ is None:
    raise ValueError("Failed to build component thermodb.")


# SECTION: build model source
# NOTE: with partially matched rules
model_source_: ComponentModelSource = build_component_model_source(
    component_thermodb=thermodb_component_,
    rules=None,
)

# =======================================
# #️⃣ THERMODB LINK CONFIGURATION
# =======================================
# NOTE: build datasource & equationsource
datasource = model_source_.data_source
equationsource = model_source_.equation_source
# log
print(datasource)
print(equationsource)

# ------------------------------------------------
# ! THERMODYNAMIC PROPERTIES
# ------------------------------------------------
# component key
component_key = selected_component.name+"-"+selected_component.state
# vapor pressure
VaPr = equationsource[component_key]['VaPr'].cal(T=300.1)
print(VaPr)

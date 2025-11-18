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

# =======================================
# #️⃣ BUILD COMPONENT THERMODB
# =======================================
# NOTE: ignore state properties
ignore_state_props = ['MW', 'VaPr', 'Cp_IG']

# SECTION: build component thermodb
# init
component_thermodb_list: List[ComponentThermoDB] = []

# build for each component
for component in components_list:
    print(f"Building ThermoDB for component: {component.name}")
    thermodb_component_ = build_component_thermodb_from_reference(
        component_name=component.name,
        component_formula=component.formula,
        component_state=component.state,
        reference_content=REFERENCE_CONTENT,
        ignore_state_props=ignore_state_props,
        component_key="Name-State"
    )

    # check
    if thermodb_component_ is None:
        logger.error(f"{component.name} ThermoDB component build failed.")
        continue

    # append to list
    component_thermodb_list.append(thermodb_component_)

# log
logger.info(f"Built {len(component_thermodb_list)} component ThermoDBs.")


# SECTION: build model source
# NOTE: with partially matched
model_source_list: List[ComponentModelSource] = build_components_model_source(
    components_thermodb=component_thermodb_list
)

# NOTE: build model source
model_source_: ModelSource = build_model_source(
    source=model_source_list
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

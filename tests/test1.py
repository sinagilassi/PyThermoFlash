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

# SECTION: bubble-temperature point calculation
# inputs
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'pressure': [101.3, 'kPa']
}

# ! bubble-point temperature calculation
res_bp = vle.bubble_temperature(inputs=inputs)
print(res_bp)

# SECTION: dew-temperature point calculation
# inputs
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'pressure': [101.3, 'kPa']
}

# ! dew-point temperature calculation
res_dp = vle.dew_temperature(inputs=inputs)
print(res_dp)

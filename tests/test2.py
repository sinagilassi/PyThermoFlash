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
# !Example 10-4, Page 452, Fundamental of Chemical Engineering Thermodynamics, Kevin D. Dahm

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: thermodb directory
thermodb_dir = os.path.join(os.getcwd(), 'tests', 'thermodb')

# ! water (H2O)
H2O_file = os.path.join(thermodb_dir, 'water-1.pkl')
H2O = ptdb.load_thermodb(H2O_file)
# check
print(H2O.check())

# general-data
res_ = H2O.select('general-data')
print(type(res_))
print(res_.table_columns)
print(res_.prop_data)

# ! ethanol (C2H5OH)
C2H5OH_file = os.path.join(thermodb_dir, 'ethanol-1.pkl')
C2H5OH = ptdb.load_thermodb(C2H5OH_file)
# check
print(C2H5OH.check())

# general-data
res_ = C2H5OH.select('general-data')
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
thub1.add_thermodb('water', H2O)
thub1.add_thermodb('ethanol', C2H5OH)

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
components = ['water', 'ethanol']

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

# SECTION: flash calculation
# inputs
inputs = {
    'mole_fraction': {'water': 0.50, 'ethanol': 0.50},
    'temperature': [30.0, 'C'],
    'pressure': [7.0, 'kPa']
}

# ! flash calculation (minimize)
res_bp = vle.flash_isothermal(inputs=inputs, solver_method='minimize')
print(res_bp)

# ! flash calculation (least_squares)
res_bp = vle.flash_isothermal(inputs=inputs)
print(res_bp)

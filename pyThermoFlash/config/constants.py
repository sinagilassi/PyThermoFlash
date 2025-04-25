# CONSTANTS
# --------------
# import libs
import math

# NOTE: eos models
PENG_ROBINSON = "PR"
REDLICH_KWONG_SOAVE = "RSK"
REDLICH_KWONG = "RK"
VAN_DER_WAALS = "VDW"

# NOTE: assumptions
RAOULT_MODEL = 'raoult'
MODIFIED_RAOULT_MODEL = 'modified-raoult'

# NOTE: activity coefficient model
VAN_LAAR_ACTIVITY_MODEL = 'van-laar'
WILSON_ACTIVITY_MODEL = 'wilson'

# NOTE: universal gas constant [J/mol.K]
R_CONST = 8.314472

# epsilon
EPS_CONST = 1e-30

# NOTE: pi
PI_CONST = math.pi

# NOTE: STP condition
# pressure [Pa]
PRESSURE_STP = 101325
# temperature [K]
TEMPERATURE_STP = 273.15
# reference temperature [K]
TEMPERATURE_REFERENCE = TEMPERATURE_STP + 25.00

# SECTION: PyThermoDBLink/PyThermoDB
DATASOURCE = "datasource"
EQUATIONSOURCE = "equationsource"

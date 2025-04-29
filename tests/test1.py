# import libs
import pyThermoFlash as ptf
from rich import print

# version
print(ptf.__version__)

# SECTION: vle model
# components
components = ['water', 'ethanol']
# init vle
vle = ptf.vle(components)
print(vle)


# SECTION: bubble point calculation
# inputs
inputs = {
    'mole_fraction': {'water': 0.5, 'ethanol': 0.5},
    'temperature': [298.15, 'K']
}

# bubble point
res_bp = vle.bubble_pressure(inputs=inputs)
print(res_bp)

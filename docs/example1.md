# 🌡️ Bubble and Dew Point Temperature

This example demonstrates the calculation of bubble and dew point temperatures for a binary mixture of benzene and toluene.

## 📋 Prerequisites

Ensure the following libraries are installed:

- `pyThermoFlash`
- `pyThermoDB`
- `pyThermoLinkDB`
- `rich`

## 🛠️ Steps

### 1️⃣ Load Thermodynamic Data

Thermodynamic data for benzene and toluene are loaded from `.pkl` files located in the `tests/thermodb` directory.

```python
# Load thermodynamic data
thermodb_dir = os.path.join(os.getcwd(), 'tests', 'thermodb')

# Benzene
C6H6_file = os.path.join(thermodb_dir, 'benzene-1.pkl')
C6H6 = ptdb.load_thermodb(C6H6_file)
print(C6H6.check())

# Toluene
C7H8_file = os.path.join(thermodb_dir, 'toluene-1.pkl')
C7H8 = ptdb.load_thermodb(C7H8_file)
print(C7H8.check())
```

### 2. Configure ThermoDB Link

Initialize the ThermoDB hub and add the components. Configure the rules using a YAML file.

```python
# Initialize ThermoDB hub
thub1 = ptdblink.init()
thub1.add_thermodb('benzene', C6H6)
thub1.add_thermodb('toluene', C7H8)

# Configure rules
thermodb_config_file = os.path.join(os.getcwd(), 'tests', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build data sources
datasource, equationsource = thub1.build()
```

### 3. Initialize VLE Model

Set up the VLE model with the components and data sources.

```python
# Components
components = ['benzene', 'toluene']

# Model source
model_source = {
    'datasource': datasource,
    'equationsource': equationsource
}

# Initialize VLE
vle = ptf.vle(components, model_source=model_source)
```

### 4. Perform Calculations

#### Bubble-Point Temperature

Calculate the bubble-point temperature using Raoult's law and the modified Raoult's law with the NRTL activity model.

```python
# Inputs
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'pressure': [101.3, 'kPa'],
}

# Activity model parameters
activity_inputs = {
    'alpha': [[0, 0.3], [0.3, 0]],
    'a_ij': [[0.0, -2.885], [2.191, 0.0]],
    'b_ij': [[0.0, 1124], [-863.7, 0.0]],
    'c_ij': [[0.0, 0.0], [0.0, 0.0]],
    'd_ij': [[0.0, 0.0], [0.0, 0.0]]
}

# Raoult's law
res_bp = vle.bubble_temperature(inputs=inputs, equilibrium_model='raoult')
print(res_bp)

# Modified Raoult's law
res_bp = vle.bubble_temperature(
    inputs=inputs,
    equilibrium_model='modified-raoult',
    activity_model='NRTL',
    activity_inputs=activity_inputs
)
print(res_bp)
```

#### Dew-Point Temperature

Calculate the dew-point temperature using Raoult's law and the modified Raoult's law with the NRTL activity model.

```python
# Raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='raoult',
    solver_method='fsolve'
)
print(res_dp)

# Modified Raoult's law
res_dp = vle.dew_temperature(
    inputs=inputs,
    equilibrium_model='modified-raoult',
    activity_model='NRTL',
    activity_inputs=activity_inputs,
    solver_method='least-squares'
)
print(res_dp)
```
# ⚡ Flash Calculation

This example demonstrates how to perform a flash calculation.

## 🛠️ Steps

### 🔹 1. Load PyThermoDB

First, load the thermodynamic data for the components involved in the calculation. In this example, we use water (H2O) and ethanol (C2H5OH):

```python
thermodb_dir = os.path.join(os.getcwd(), 'tests', 'thermodb')

# Load water data
H2O_file = os.path.join(thermodb_dir, 'water-1.pkl')
H2O = ptdb.load_thermodb(H2O_file)
print(H2O.check())

# Load ethanol data
C2H5OH_file = os.path.join(thermodb_dir, 'ethanol-1.pkl')
C2H5OH = ptdb.load_thermodb(C2H5OH_file)
print(C2H5OH.check())
```

### 🔹 2. Configure PyThermoDBLink

Initialize the ThermoDB hub and configure the rules for linking the thermodynamic data:

```python
thub1 = ptdblink.init()
thub1.add_thermodb('water', H2O)
thub1.add_thermodb('ethanol', C2H5OH)

thermodb_config_file = os.path.join(os.getcwd(), 'tests', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

datasource, equationsource = thub1.build()
```

### 🔹 3. Initialize VLE Model

Set up the vapor-liquid equilibrium (VLE) model with the required components and sources:

```python
components = ['water', 'ethanol']
model_source = {
    'datasource': datasource,
    'equationsource': equationsource
}

vle = ptf.vle(components, model_source=model_source)
```

### 🔹 4. Perform Flash Calculation

Define the input conditions and perform the flash calculation:

```python
inputs = {
    'mole_fraction': {'water': 0.50, 'ethanol': 0.50},
    'temperature': [30.0, 'C'],
    'pressure': [7.0, 'kPa']
}

res_bp = vle.flash_isothermal(
    inputs=inputs,
)
print(res_bp)
```

This will output the results of the flash calculation, including the equilibrium state of the system.
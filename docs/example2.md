# 🫧 Bubble and Dew Point Pressure

This example demonstrates how to perform bubble and dew point pressure calculations.

## 📂 Loading Thermodynamic Data

Thermodynamic data for benzene and toluene are loaded from `.pkl` files located in the `tests/thermodb` directory:

```python
thermodb_dir = os.path.join(os.getcwd(), 'tests', 'thermodb')
C6H6_file = os.path.join(thermodb_dir, 'benzene-1.pkl')
C6H6 = ptdb.load_thermodb(C6H6_file)
C7H8_file = os.path.join(thermodb_dir, 'toluene-1.pkl')
C7H8 = ptdb.load_thermodb(C7H8_file)
```

## 🔗 Configuring ThermoDB Link

A `thermodb` hub is initialized, and the components are added to it. A configuration file is used to define rules for the `thermodb`:

```python
thub1 = ptdblink.init()
thub1.add_thermodb('benzene', C6H6)
thub1.add_thermodb('toluene', C7H8)
thermodb_config_file = os.path.join(os.getcwd(), 'tests', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)
```

## ⚙️ Initializing the VLE Model

The vapor-liquid equilibrium (VLE) model is initialized with the components and the data sources:

```python
components = ['benzene', 'toluene']
model_source = {
    'datasource': datasource,
    'equationsource': equationsource
}
vle = ptf.vle(components, model_source=model_source)
```

## Bubble Point Pressure Calculation

Bubble point pressure is calculated using Raoult's law and a modified Raoult's law with the NRTL activity model:

```python
inputs = {
    'mole_fraction': {'benzene': 0.26, 'toluene': 0.74},
    'temperature': [80, 'C'],
}

# Raoult's law
res_bp = vle.bubble_pressure(inputs=inputs, equilibrium_model='raoult')
print(res_bp)

# Modified Raoult's law
activity_inputs = {
    'alpha': [[0, 0.3], [0.3, 0]],
    'a_ij': [[0.0, -2.885], [2.191, 0.0]],
    'b_ij': [[0.0, 1124], [-863.7, 0.0]],
    'c_ij': [[0.0, 0.0], [0.0, 0.0]],
    'd_ij': [[0.0, 0.0], [0.0, 0.0]]
}
res_bp = vle.bubble_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_inputs=activity_inputs
)
print(res_bp)
```

## Dew Point Pressure Calculation

Dew point pressure is calculated similarly:

```python
# Raoult's law
res_dp = vle.dew_pressure(inputs=inputs, equilibrium_model='raoult')
print(res_dp)

# Modified Raoult's law
res_dp = vle.dew_pressure(
    inputs=inputs, equilibrium_model='modified-raoult',
    activity_model='NRTL', activity_inputs=activity_inputs
)
print(res_dp)
```
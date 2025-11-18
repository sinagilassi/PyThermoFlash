# import libs
from typing import Generic, TypeVar, Annotated, List, Optional, Dict, Any
from pydantic import BaseModel, Field, ConfigDict, BeforeValidator
from pythermodb_settings.models import Pressure, Temperature, Component

T = TypeVar('T')


class Quantity(BaseModel, Generic[T]):
    """Generic quantity with value + unit (e.g. Pa, K, s, dimensionless)."""
    value: Optional[T]
    unit: str
    symbol: str


class ComponentProps(BaseModel, Generic[T]):
    component: Component
    properties: Dict[str, Quantity[T]]


# SECTION: pressure result model

# {
#     'bubble_pressure': {'value': 54880.72088600001, 'unit': 'Pa'},
#     'temperature': {'value': 353.15, 'unit': 'K'},
#     'feed_mole_fraction': array([0.26, 0.74]),
#     'vapor_mole_fraction': array([0.47787172, 0.52212828]),
#     'liquid_mole_fraction': array([0.26, 0.74]),
#     'mole_fraction_sum': {'zi': 1.0, 'xi': 1.0, 'yi': 0.9999999999999999},
#     'vapor_pressure': {'value': array([100869.017 ,  38722.6709]), 'unit': 'Pa'},
#     'activity_coefficient': {'value': array([1., 1.]), 'unit': 'dimensionless'},
#     'K_ratio': {'value': array([1.83796815, 0.70557876]), 'unit': 'dimensionless'},
#     'message': 'Bubble Pressure Calculation',
#     'components': ['benzene-l', 'toluene-l'],
#     'equilibrium_model': 'raoult',
#     'fugacity_model': None,
#     'activity_model': None,
#     'computation_time': {'value': 0.0009999275207519531, 'unit': 's'}
# }

# {
#     'dew_pressure': {'value': 46108.761202690344, 'unit': 'Pa'},
#     'temperature': {'value': {'value': 353.15, 'unit': 'K'}, 'unit': 'K'},
#     'feed_mole_fraction': array([0.26, 0.74]),
#     'vapor_mole_fraction': array([0.26, 0.74]),
#     'liquid_mole_fraction': array([0.11884995, 0.88115005]),
#     'mole_fraction_sum': {'zi': 1.0, 'xi': 1.0, 'yi': 1.0},
#     'vapor_pressure': {'value': array([100869.017 ,  38722.6709]), 'unit': 'Pa'},
#     'activity_coefficient': {'value': array([1., 1.]), 'unit': 'dimensionless'},
#     'K_ratio': {'value': array([2.18763234, 0.83981156]), 'unit': 'dimensionless'},
#     'max_iter': 500,
#     'iteration': 0,
#     'tolerance': 1e-06,
#     'message': 'Dew Pressure Calculation',
#     'components': ['benzene-l', 'toluene-l'],
#     'equilibrium_model': 'raoult',
#     'fugacity_model': None,
#     'activity_model': None,
#     'computation_time': {'value': 0.0010080337524414062, 'unit': 's'}
# }


class BubblePressureResult(BaseModel):
    model_config = ConfigDict(
        extra="allow"
    )

    bubble_pressure: Pressure = Field(
        ...,
        description="Calculated bubble-point pressure."
    )
    temperature: Temperature = Field(
        ...,
        description="System temperature."
    )
    component_ids: List[str] = Field(
        ...,
        description="List of component IDs involved in the calculation."
    )

    feed_mole_fraction: List[Component]
    vapor_mole_fraction: List[Component]
    liquid_mole_fraction: List[Component]
    mole_fraction_sum: Dict[str, float]

    component_props: List[ComponentProps] = Field(
        ...,
        description="List of component properties used in the calculation."
    )

    message: str

    equilibrium_model: str | None
    fugacity_model: str | None
    activity_model: str | None

    computation_time: Quantity[float] | None


class DewPressureResult(BaseModel):
    model_config = ConfigDict(
        extra="allow"
    )

    dew_pressure: Pressure = Field(
        ...,
        description="Calculated dew-point pressure."
    )
    temperature: Temperature = Field(
        ...,
        description="System temperature."
    )
    component_ids: List[str] = Field(
        ...,
        description="List of component IDs involved in the calculation."
    )

    feed_mole_fraction: List[Component]
    vapor_mole_fraction: List[Component]
    liquid_mole_fraction: List[Component]
    mole_fraction_sum: Dict[str, float]

    component_props: List[ComponentProps] = Field(
        ...,
        description="List of component properties used in the calculation."
    )

    message: str

    equilibrium_model: str | None
    fugacity_model: str | None
    activity_model: str | None

    computation_time: Quantity[float] | None

    max_iter: int | None
    iteration: int | None
    tolerance: float | None

# SECTION: temperature result model

# {
#     'bubble_temperature': {'value': 373.07718645402184, 'unit': 'K'},
#     'pressure': {'value': 101299.99999999999, 'unit': 'Pa'},
#     'feed_mole_fraction': [0.26, 0.74],
#     'liquid_mole_fraction': [0.26, 0.74],
#     'vapor_mole_fraction': [0.4605291602566635, 0.53947083994077],
#     'mole_fraction_sum': {'zi': 1.0, 'xi': 1.0, 'yi': 1.0000000001974336},
#     'vapor_pressure': {'value': [179429.2459, 73849.1839], 'unit': 'Pa'},
#     'K_ratio': {'value': [1.7712660009871672, 0.7290146485686082], 'unit': 'dimensionless'},
#     'activity_coefficient': {'value': [1.0, 1.0], 'unit': 'dimensionless'},
#     'message': 'Bubble Temperature Calculation',
#     'components': ['benzene-l', 'toluene-l'],
#     'equilibrium_model': 'raoult',
#     'fugacity_model': None,
#     'activity_model': None,
#     'solver_method': 'root',
#     'computation_time': {'value': 0.005999326705932617, 'unit': 's'}
# }

# {
#     'dew_temperature': {'value': 378.1633323578604, 'unit': 'K'},
#     'pressure': {'value': 101299.99999999999, 'unit': 'Pa'},
#     'feed_mole_fraction': [0.26, 0.74],
#     'liquid_mole_fraction': [0.1281658714162829, 0.8718341283338952],
#     'vapor_mole_fraction': [0.26, 0.74],
#     'mole_fraction_sum': {'xi': 0.9999999997501781, 'yi': 1.0, 'zi': 1.0},
#     'vapor_pressure': {'value': [205499.3245, 85981.9518], 'unit': 'Pa'},
#     'K_ratio': {'value': [0.4929456592933958, 1.1781542274782366], 'unit': 'dimensionless'},
#     'activity_coefficient': {'value': [1.0, 1.0], 'unit': 'dimensionless'},
#     'message': 'Dew Temperature Calculation',
#     'components': ['benzene-l', 'toluene-l'],
#     'equilibrium_model': 'raoult',
#     'fugacity_model': None,
#     'activity_model': None,
#     'solver_method': 'root',
#     'computation_time': {'value': 0.012703895568847656, 'unit': 's'}
# }


class BubbleTemperatureResult(BaseModel):
    model_config = ConfigDict(
        extra="allow"
    )

    bubble_temperature: Temperature = Field(
        ...,
        description="Calculated bubble-point temperature."
    )
    pressure: Pressure = Field(
        ...,
        description="System pressure."
    )
    component_ids: List[str] = Field(
        ...,
        description="List of component IDs involved in the calculation."
    )

    feed_mole_fraction: List[Component]
    vapor_mole_fraction: List[Component]
    liquid_mole_fraction: List[Component]
    mole_fraction_sum: Dict[str, float]

    component_props: List[ComponentProps] = Field(
        ...,
        description="List of component properties used in the calculation."
    )

    message: str

    equilibrium_model: str | None
    fugacity_model: str | None
    activity_model: str | None

    solver_method: str | None

    computation_time: Quantity[float] | None


class DewTemperatureResult(BaseModel):
    model_config = ConfigDict(
        extra="allow"
    )

    dew_temperature: Temperature = Field(
        ...,
        description="Calculated dew-point temperature."
    )
    pressure: Pressure = Field(
        ...,
        description="System pressure."
    )
    component_ids: List[str] = Field(
        ...,
        description="List of component IDs involved in the calculation."
    )

    feed_mole_fraction: List[Component]
    vapor_mole_fraction: List[Component]
    liquid_mole_fraction: List[Component]
    mole_fraction_sum: Dict[str, float]

    component_props: List[ComponentProps] = Field(
        ...,
        description="List of component properties used in the calculation."
    )

    message: str

    equilibrium_model: str | None
    fugacity_model: str | None
    activity_model: str | None

    solver_method: str | None

    computation_time: Quantity[float] | None

# SECTION: flash result model

# {
#     'message': 'Is Flashable Calculation',
#     'pressure': {'value': 7000.000000000001, 'unit': 'Pa'},
#     'temperature': {'value': 303.15, 'unit': 'K'},
#     'feed_mole_fraction': {'water-l': 0.5, 'ethanol-l': 0.5},
#     'components': ['water-l', 'ethanol-l'],
#     'equilibrium_model': 'raoult',
#     'flash_checker_res': True,
#     'computation_time': {'value': 0.0022618770599365234, 'unit': 's'}
# }


class FlashIsothermalResult(BaseModel):
    model_config = ConfigDict(
        extra="allow"
    )

    message: str

    pressure: Pressure = Field(
        ...,
        description="System pressure."
    )
    temperature: Temperature = Field(
        ...,
        description="System temperature."
    )
    component_ids: List[str] = Field(
        ...,
        description="List of component IDs involved in the calculation."
    )

    feed_mole_fraction: List[Component]

    equilibrium_model: str | None

    flash_checker_res: bool | None

    computation_time: Quantity[float] | None

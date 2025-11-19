# import libs
import logging
from typing import Dict, Any, List, Union, Literal
from pythermodb_settings.models import Pressure, Temperature, Component
# local
from ..models import (
    BubblePressureResult,
    DewPressureResult,
    BubbleTemperatureResult,
    DewTemperatureResult,
    ComponentProps,
    Quantity,
    FlashIsothermalResult,
    CheckFlashIsothermalResult
)

# NOTE: setup logger
logger = logging.getLogger(__name__)


def _prepare_pressure_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
        model_type: Literal['BP', 'DP', 'BT', 'DT', 'IF', 'CIF'],
) -> Union[
    BubblePressureResult,
    DewPressureResult,
    BubbleTemperatureResult,
    DewTemperatureResult,
    FlashIsothermalResult,
    CheckFlashIsothermalResult,
    None
]:
    '''
    Prepare result structure for output.
    '''
    try:
        # SECTION: parse result
        # NOTE: temperature
        temperature_ = data.get('temperature', None)
        if temperature_ is not None:
            temperature = Temperature(
                value=float(temperature_['value']),
                unit=temperature_['unit']
            )
        else:
            temperature = Temperature(
                value=0.0,
                unit='K'
            )

        # NOTE: pressure
        pressure_ = data.get('pressure', None)
        if pressure_ is not None:
            pressure = Pressure(
                value=float(pressure_['value']),
                unit=pressure_['unit']
            )
        else:
            pressure = Pressure(
                value=0.0,
                unit='Pa'
            )

        # NOTE: pressure calc result
        calc_pressure = data.get(
            'bubble_pressure', None) or data.get('dew_pressure', None)
        if calc_pressure is not None:
            result_pressure = Pressure(
                value=float(calc_pressure['value']),
                unit=calc_pressure['unit']
            )
        else:
            result_pressure = Pressure(
                value=0.0,
                unit='Pa'
            )

        # NOTE: temperature calc result
        calc_temperature = data.get(
            'bubble_temperature', None) or data.get('dew_temperature', None)

        if calc_temperature is not None:
            result_temperature = Temperature(
                value=float(calc_temperature['value']),
                unit=calc_temperature['unit']
            )
        else:
            result_temperature = Temperature(
                value=0.0,
                unit='K'
            )

        # NOTE:component ids
        component_ids = data.get('components', [])

        # NOTE: mole fractions
        feed_mole_fraction_ = data.get('feed_mole_fraction', [])
        vapor_mole_fraction_ = data.get('vapor_mole_fraction', [])
        liquid_mole_fraction_ = data.get('liquid_mole_fraction', [])
        mole_fraction_sum = data.get('mole_fraction_sum', [])

        # NOTE: activity coefficients, K ratios, vapor pressures
        activity_coefficients = data.get('activity_coefficient', None)
        K_ratios = data.get('K_ratio', None)
        vapor_pressures = data.get('vapor_pressure', None)

        # NOTE: message
        message = data.get('message', 'Not specified.')

        # NOTE: models
        equilibrium_model = data.get('equilibrium_model', None)
        fugacity_model = data.get('fugacity_model', None)
        activity_model = data.get('activity_model', None)

        # NOTE: computation time
        computation_time = data.get('computation_time', None)
        # >> check
        if computation_time is not None:
            computation_time = Quantity[float](
                value=float(computation_time['value']),
                unit=computation_time['unit'],
                symbol='t'
            )

        # NOTE: solver method
        solver_method = data.get('solver_method', None)

        # SECTION: feed
        feed_mole_fraction: List[Component] = []
        vapor_mole_fraction: List[Component] = []
        liquid_mole_fraction: List[Component] = []

        # iterate components
        for i, comp_name in enumerate(component_names):
            # feed mole fraction
            mole_frac_value = feed_mole_fraction_[i] if i < len(
                feed_mole_fraction_) else None
            # vapor mole fraction
            vapor_mole_frac_value = vapor_mole_fraction_[i] if i < len(
                vapor_mole_fraction_) else None
            # liquid mole fraction
            liquid_mole_frac_value = liquid_mole_fraction_[i] if i < len(
                liquid_mole_fraction_) else None

            comp = next(
                (c for c in components if c.name == comp_name),
                None
            )

            # NOTE: feed mole fraction
            if (
                comp is not None and
                mole_frac_value is not None
            ):
                # >> append
                feed_mole_fraction.append(
                    Component(
                        name=comp.name,
                        formula=comp.formula,
                        state=comp.state,
                        mole_fraction=float(mole_frac_value)
                    )
                )

            # NOTE: vapor mole fraction
            if (
                comp is not None and
                vapor_mole_frac_value is not None
            ):
                # >> append
                vapor_mole_fraction.append(
                    Component(
                        name=comp.name,
                        formula=comp.formula,
                        state=comp.state,
                        mole_fraction=float(vapor_mole_frac_value)
                    )
                )

            # NOTE: liquid mole fraction
            if (
                comp is not None and
                liquid_mole_frac_value is not None
            ):
                # >> append
                liquid_mole_fraction.append(
                    Component(
                        name=comp.name,
                        formula=comp.formula,
                        state=comp.state,
                        mole_fraction=float(liquid_mole_frac_value)
                    )
                )

        # NOTE: calculation setup
        max_iter = data.get('max_iter', None)
        tolerance = data.get('tolerance', None)
        iteration = data.get('iteration', None)

        # NOTE: flash properties
        V_F_ratio = data.get('V_F_ratio', None)
        # >> check
        if V_F_ratio is not None:
            V_F_ratio = Quantity[float](
                value=float(V_F_ratio['value']),
                unit=V_F_ratio['unit'],
                symbol=V_F_ratio['symbol']
            )
        else:
            V_F_ratio = Quantity[float](
                value=0.0,
                unit='-',
                symbol='V/F'
            )

        L_F_ratio = data.get('L_F_ratio', None)
        # >> check
        if L_F_ratio is not None:
            L_F_ratio = Quantity[float](
                value=float(L_F_ratio['value']),
                unit=L_F_ratio['unit'],
                symbol=L_F_ratio['symbol']
            )
        else:
            L_F_ratio = Quantity[float](
                value=0.0,
                unit='-',
                symbol='L/F'
            )

        solver_message = data.get('solver_message', None)
        flash_checker = data.get('flash_checker', None)
        flash_checker_res = data.get('flash_checker_res', None)

        # SECTION: component properties
        component_props: List[ComponentProps] = []

        for i, comp_name in enumerate(component_names):
            # find component
            component = next(
                (comp for comp in components if comp.name == comp_name),
                None
            )

            # >> check
            if component is None:
                logger.warning(
                    f"Component '{comp_name}' not found in components list.")
                continue

            comp_prop = ComponentProps(
                component=component,
                properties={}
            )

            # >> vapor pressure
            if vapor_pressures is not None:
                vp_value = vapor_pressures['value'][i]
                vp_unit = vapor_pressures['unit']
                comp_prop.properties['vapor_pressure'] = Quantity[float](
                    value=float(vp_value),
                    unit=vp_unit,
                    symbol='VaPr'
                )

            # >> activity coefficient
            if activity_coefficients is not None:
                ac_value = activity_coefficients['value'][i]
                ac_unit = activity_coefficients['unit']
                comp_prop.properties['activity_coefficient'] = Quantity[float](
                    value=float(ac_value),
                    unit=ac_unit,
                    symbol='AcCo'
                )

            # >> K ratio
            if K_ratios is not None:
                K_value = K_ratios['value'][i]
                K_unit = K_ratios['unit']
                comp_prop.properties['K_ratio'] = Quantity[float](
                    value=float(K_value),
                    unit=K_unit,
                    symbol='Kxy'
                )

            # NOTE: append other properties as needed
            component_props.append(comp_prop)

        # SECTION: init result
        if model_type == 'BP':
            # ! Bubble Pressure Result
            result_structure = BubblePressureResult(
                bubble_pressure=result_pressure,
                temperature=temperature,
                component_ids=component_ids,
                feed_mole_fraction=feed_mole_fraction,
                vapor_mole_fraction=vapor_mole_fraction,
                liquid_mole_fraction=liquid_mole_fraction,
                mole_fraction_sum=mole_fraction_sum,
                component_props=component_props,
                message=message,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                computation_time=computation_time
            )
        elif model_type == 'DP':
            # ! Dew Pressure Result
            result_structure = DewPressureResult(
                dew_pressure=result_pressure,
                temperature=temperature,
                component_ids=component_ids,
                feed_mole_fraction=feed_mole_fraction,
                vapor_mole_fraction=vapor_mole_fraction,
                liquid_mole_fraction=liquid_mole_fraction,
                mole_fraction_sum=mole_fraction_sum,
                component_props=component_props,
                message=message,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                computation_time=computation_time,
                max_iter=max_iter,
                tolerance=tolerance,
                iteration=iteration
            )
        elif model_type == 'BT':
            # ! Bubble Temperature Result
            result_structure = BubbleTemperatureResult(
                bubble_temperature=result_temperature,
                pressure=pressure,
                component_ids=component_ids,
                feed_mole_fraction=feed_mole_fraction,
                vapor_mole_fraction=vapor_mole_fraction,
                liquid_mole_fraction=liquid_mole_fraction,
                mole_fraction_sum=mole_fraction_sum,
                component_props=component_props,
                message=message,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                solver_method=solver_method,
                computation_time=computation_time
            )
        elif model_type == 'DT':
            # ! Dew Temperature Result
            result_structure = DewTemperatureResult(
                dew_temperature=result_temperature,
                pressure=pressure,
                component_ids=component_ids,
                feed_mole_fraction=feed_mole_fraction,
                vapor_mole_fraction=vapor_mole_fraction,
                liquid_mole_fraction=liquid_mole_fraction,
                mole_fraction_sum=mole_fraction_sum,
                component_props=component_props,
                message=message,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                solver_method=solver_method,
                computation_time=computation_time
            )
        elif model_type == 'IF':
            # ! Flash Isothermal Result
            result_structure = FlashIsothermalResult(
                V_F_ratio=V_F_ratio,
                L_F_ratio=L_F_ratio,
                temperature=temperature,
                pressure=pressure,
                component_ids=component_ids,
                feed_mole_fraction=feed_mole_fraction,
                vapor_mole_fraction=vapor_mole_fraction,
                liquid_mole_fraction=liquid_mole_fraction,
                mole_fraction_sum=mole_fraction_sum,
                component_props=component_props,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                flash_checker=flash_checker,
                flash_checker_res=flash_checker_res,
                solver_method=solver_method,
                solver_message=solver_message,
                message=message,
                computation_time=computation_time,
            )
        elif model_type == 'CIF':
            # ! Check Flash Isothermal Result
            result_structure = CheckFlashIsothermalResult(
                pressure=pressure,
                temperature=temperature,
                component_ids=component_ids,
                components=components,
                feed_mole_fraction=feed_mole_fraction,
                equilibrium_model=equilibrium_model,
                flash_checker_res=flash_checker_res,
                message=message,
                computation_time=computation_time,
            )
        else:
            result_structure = None

        return result_structure
    except Exception as e:
        logger.error(f"Error preparing result structure: {e}")
        raise e


def prepare_bubble_pressure_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> BubblePressureResult | None:
    '''
    Prepare bubble pressure result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    BubblePressureResult | None
        Structured BubblePressureResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='BP',
    )

    if isinstance(res, BubblePressureResult):
        return res
    return None


def prepare_dew_pressure_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> DewPressureResult | None:
    '''
    Prepare dew pressure result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    DewPressureResult | None
        Structured DewPressureResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='DP',
    )

    if isinstance(res, DewPressureResult):
        return res
    return None


def prepare_bubble_temperature_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> BubbleTemperatureResult | None:
    '''
    Prepare bubble temperature result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    BubbleTemperatureResult | None
        Structured BubbleTemperatureResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='BT',
    )

    if isinstance(res, BubbleTemperatureResult):
        return res
    return None


def prepare_dew_temperature_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> DewTemperatureResult | None:
    '''
    Prepare dew temperature result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    DewTemperatureResult | None
        Structured DewTemperatureResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='DT',
    )

    if isinstance(res, DewTemperatureResult):
        return res
    return None


def prepare_flash_isothermal_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> FlashIsothermalResult | None:
    '''
    Prepare flash isothermal result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    FlashIsothermalResult | None
        Structured FlashIsothermalResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='IF',
    )

    if isinstance(res, FlashIsothermalResult):
        return res
    return None


def prepare_check_flash_isothermal_result_structure(
        component_names: List[str],
        components: List[Component],
        data: dict,
) -> CheckFlashIsothermalResult | None:
    '''
    Prepare check flash isothermal result structure for output.

    Parameters
    ----------
    component_names : List[str]
        List of component names involved in the calculation.
    components : List[Component]
        List of Component objects used in the calculation.
    data : dict
        Raw data dictionary containing calculation results.

    Returns
    -------
    CheckFlashIsothermalResult | None
        Structured CheckFlashIsothermalResult object or None if preparation fails.
    '''
    res = _prepare_pressure_result_structure(
        component_names=component_names,
        components=components,
        data=data,
        model_type='CIF',
    )

    if isinstance(res, CheckFlashIsothermalResult):
        return res
    return None

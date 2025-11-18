# import libs
import logging
from typing import Any, Dict, Literal, List, Optional
from pythermodb_settings.models import Component, Temperature, Pressure
from pythermodb_settings.utils import set_component_id
from pyThermoLinkDB.models import ModelSource
# locals
from ..utils import set_feed_specification
from ..docs.vle import VLE

# NOTE: logger
logger = logging.getLogger(__name__)

# SECTION: Input generator


def _input_generator(
    model_type: Literal[
        'bpp',
        'dpp',
        'bpt',
        'dpt',
        'if'
    ],
    components: List[Component],
    model_source: ModelSource,
    temperature: Optional[Temperature] = None,
    pressure: Optional[Pressure] = None,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
):
    '''
    Input generator for thermodynamic property calculations.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'vle_model': Initialized VLE model.
        - 'inputs': Dictionary of inputs for the calculation.

    Notes
    -----
    - Validates inputs for components, temperature, pressure, and model source.
    - Prepares component IDs and feed specifications based on provided components.
    - Initializes a VLE model with the prepared inputs.
    - model_type: 'bpp' for bubble-point pressure, 'dpp' for dew-point pressure,
    'bpt' for bubble-point temperature, 'dpt' for dew-point temperature, 'if' for isothermal flash.
    '''
    try:
        # SECTION: validate inputs
        # ! component
        if not isinstance(components, list):
            raise ValueError(
                "Invalid component input. Must be a Component object.")

        # >> validate each component
        for comp in components:
            if not isinstance(comp, Component):
                raise ValueError(
                    "Invalid component input. Must be a Component object.")

        # ! temperature
        if temperature is not None:
            if not isinstance(temperature, Temperature):
                raise ValueError(
                    "Invalid temperature input. Must be a Temperature object.")

        # ! pressure
        if pressure is not None:
            if not isinstance(pressure, Pressure):
                raise ValueError(
                    "Invalid pressure input. Must be a Pressure object.")

        # ! model source
        if not isinstance(model_source, ModelSource):
            raise ValueError(
                "Invalid model_source input. Must be a ModelSource object.")

        # SECTION: input preparation
        try:
            # NOTE component id configuration
            # set component id
            component_ids: List[str] = [
                set_component_id(
                    component=comp,
                    component_key=component_key
                ) for comp in components
            ]

            # NOTE: set feed specification
            feed_spec: Dict[str, float] = set_feed_specification(
                components=components,
                component_key=component_key
            )

            # model input
            inputs: Dict[str, Any] = {
                'mole_fraction': feed_spec,
            }

            # >> temperature input
            if model_type == 'bpp':
                if temperature is None:
                    raise ValueError(
                        "Temperature input is required for bubble-point pressure calculation.")
                inputs["temperature"] = [temperature.value, temperature.unit]

            # >> pressure input
            if model_type == 'bpt':
                if pressure is None:
                    raise ValueError(
                        "Pressure input is required for bubble-point temperature calculation.")
                inputs["pressure"] = [pressure.value, pressure.unit]

            # model source
            model_source_dict = {
                "datasource": model_source.data_source,
                "equationsource": model_source.equation_source
            }
        except Exception as e:
            logger.error(f"Input preparation failed!, {e}")
            raise

        # SECTION: init VLE model
        try:
            vle_model = VLE(
                components=component_ids,
                model_source=model_source_dict,
            )
        except Exception as e:
            logger.error(f"VLE model initialization failed!, {e}")
            raise

        # SECTION: return inputs
        return {
            'vle_model': vle_model,
            'inputs': inputs
        }
    except Exception as e:
        logger.error(f"Error in _pc_input_generator: {e}")
        raise e

# SECTION: Bubble-point pressure (BPP) calculation


def calc_bubble_point_pressure(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
    equilibrium_model: Literal[
        'raoult',
        'modified-raoult'
    ] = 'raoult',
    fugacity_model: Optional[Literal[
        'vdW', 'PR', 'RK', 'SRK'
    ]] = None,
    activity_model: Optional[Literal[
        'NRTL', 'UNIQUAC'
    ]] = None,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
    message: str | None = None,
    **kwargs
):
    '''
    The Bubble-Point Pressure (BPP) calculation determines the pressure at which the first bubble of vapor forms when a liquid mixture is heated at a constant temperature. It is used to find the pressure for a given temperature at which the liquid will begin to vaporize.

    Parameters
    ----------
    components : List[Component]
        List of Component objects representing the mixture.
    temperature : Temperature
        Temperature object specifying the temperature at which to calculate the bubble-point pressure.
    model_source : ModelSource
        ModelSource object containing the data source and equation of state information.
    equilibrium_model : Literal['raoult', 'modified-raoult'], optional
        The equilibrium model to use for the calculation. Options are 'raoult' or 'modified-raoult'. Default is 'raoult'.
    fugacity_model : Optional[Literal['vdW', 'PR', 'RK', 'SRK']], optional
        The fugacity model to use for the calculation. Options are 'vdW' (Van der Waals), 'PR' (Peng-Robinson), 'RK' (Redlich-Kwong), or 'SRK' (Soave-Redlich-Kwong). Default is None.
    activity_model : Optional[Literal['NRTL', 'UNIQUAC']], optional
        The activity coefficient model to use for the calculation. Options are 'NRTL' (Non-Random Two-Liquid) or 'UNIQUAC' (Universal Quasi-Chemical). Default is None.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State', 'Formula-Name-State'], optional
        The key to identify components. Options are "Name-State" or "Formula-State". Default is "Name-State".
    message : str | None, optional
        Custom message to display during the calculation. Default is None.
    **kwargs : dict, optional
        Additional keyword arguments specific to the calculation.
        - activity_inputs : dict, model parameters.
        - activity_coefficients : dict, activity coefficients for each component.

    Returns
    -------
    res : dict
        Dictionary containing the results of the calculation.
        - bubble_pressure: bubble pressure [Pa]
        - temperature: temperature [K]
        - feed_mole_fraction: liquid mole fraction (xi)
        - vapor_mole_fraction: vapor mole fraction (yi)
        - liquid_mole_fraction: liquid mole fraction (xi)
        - vapor_pressure: vapor pressure [Pa]
        - activity_coefficient: activity coefficient [dimensionless]
        - K_ratio: K-ratio [dimensionless]
        - equilibrium_model: equilibrium model used for the calculation
        - fugacity_model: fugacity model used for the calculation
        - activity_model: activity coefficient model used for the calculation
        - components: list of components used in the calculation
        - message: message displayed during the calculation
        - computation_time: computation time [s]

    Notes
    -----
    - Temperature must be in Kelvin, Celsius, or Fahrenheit, and will be converted to Kelvin.
    - Mole fractions must be between 0 and 1.
    - The function will raise a ValueError if the inputs are not valid.
    - The default equilibrium model is 'raoult', and the default activity model is 'NRTL'.
    - The function will return the bubble pressure in Pascals.
    - The activity coefficient for each component can be directly defined in the activity_inputs dictionary.

    '''
    try:
        # SECTION: generate inputs
        input_dict = _input_generator(
            model_type='bpp',
            components=components,
            model_source=model_source,
            temperature=temperature,
            component_key=component_key,
        )

        # unpack inputs
        vle_model: VLE = input_dict['vle_model']
        inputs: Dict[str, Any] = input_dict['inputs']

        # SECTION: perform BPP calculation
        try:
            res = vle_model.bubble_pressure(
                inputs=inputs,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                message=message,
                **kwargs
            )
            return res
        except Exception as e:
            logger.error(f"BPP calculation failed!, {e}")
            raise

    except Exception as e:
        logger.error(f"Error in calc_bpp: {e}")
        raise e

# SECTION: Dew-point pressure (DPP) calculation


def calc_dew_point_pressure(
    components: List[Component],
    temperature: Temperature,
    model_source: ModelSource,
    equilibrium_model: Literal[
        'raoult',
        'modified-raoult'
    ] = 'raoult',
    fugacity_model: Optional[Literal[
        'vdW', 'PR', 'RK', 'SRK'
    ]] = None,
    activity_model: Optional[Literal[
        'NRTL', 'UNIQUAC'
    ]] = None,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
    message: str | None = None,
    **kwargs
):
    '''
    The dew-point pressure (DPP) calculation determines the pressure at which the first drop of liquid condenses when a vapor mixture is cooled at a constant temperature. It is used to find the pressure at which vapor will begin to condense.

    Parameters
    ----------
    components : List[Component]
        List of Component objects representing the mixture.
    temperature : Temperature
        Temperature object specifying the temperature at which to calculate the bubble-point pressure.
    model_source : ModelSource
        ModelSource object containing the data source and equation of state information.
    equilibrium_model : Literal['raoult', 'modified-raoult'], optional
        The equilibrium model to use for the calculation. Options are 'raoult' or 'modified-raoult'. Default is 'raoult'.
    fugacity_model : Optional[Literal['vdW', 'PR', 'RK', 'SRK']], optional
        The fugacity model to use for the calculation. Options are 'vdW' (Van der Waals), 'PR' (Peng-Robinson), 'RK' (Redlich-Kwong), or 'SRK' (Soave-Redlich-Kwong). Default is None.
    activity_model : Optional[Literal['NRTL', 'UNIQUAC']], optional
        The activity coefficient model to use for the calculation. Options are 'NRTL' (Non-Random Two-Liquid) or 'UNIQUAC' (Universal Quasi-Chemical). Default is None.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State', 'Formula-Name-State'], optional
        The key to identify components. Options are "Name-State" or "Formula-State". Default is "Name-State".
    message : str | None, optional
        Custom message to display during the calculation. Default is None.
    **kwargs : dict, optional
        Additional keyword arguments specific to the calculation.
        - activity_inputs : dict, model parameters.
        - activity_coefficients : dict, activity coefficients for each component.

    Returns
    -------
    res : dict
        Dictionary containing the results of the calculation.
        - dew_pressure: dew pressure [Pa]
        - temperature: temperature [K]
        - feed_mole_fraction: liquid mole fraction (xi)
        - vapor_mole_fraction: vapor mole fraction (yi)
        - liquid_mole_fraction: liquid mole fraction (xi)
        - vapor_pressure: vapor pressure [Pa]
        - activity_coefficient: activity coefficient [dimensionless]
        - K_ratio: K-ratio [dimensionless]
        - equilibrium_model: equilibrium model used for the calculation
        - fugacity_model: fugacity model used for the calculation
        - activity_model: activity coefficient model used for the calculation
        - components: list of components used in the calculation
        - message: message displayed during the calculation
        - computation_time: computation time [s]

    Notes
    -----
    - Temperature must be in Kelvin, Celsius, or Fahrenheit, and will be converted to Kelvin.
    - Mole fractions must be between 0 and 1.
    - The function will raise a ValueError if the inputs are not valid.
    - The default equilibrium model is 'raoult', and the default activity model is 'NRTL'.
    - The function will return the bubble pressure in Pascals.
    - The activity coefficient for each component can be directly defined in the activity_inputs dictionary.
    '''
    try:
        # SECTION: generate inputs
        input_dict = _input_generator(
            model_type='dpp',
            components=components,
            model_source=model_source,
            temperature=temperature,
            component_key=component_key,
        )

        # unpack inputs
        vle_model: VLE = input_dict['vle_model']
        inputs: Dict[str, Any] = input_dict['inputs']

        # SECTION: perform BPP calculation
        try:
            res = vle_model.dew_pressure(
                inputs=inputs,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                message=message,
                **kwargs
            )
            return res
        except Exception as e:
            logger.error(f"DPP calculation failed!, {e}")
            raise

    except Exception as e:
        logger.error(f"Error in calc_dpp: {e}")
        raise e


# SECTION: Bubble-point temperature (BPT) calculation


def calc_bubble_point_temperature(
    components: List[Component],
    pressure: Pressure,
    model_source: ModelSource,
    equilibrium_model: Literal[
        'raoult',
        'modified-raoult'
    ] = 'raoult',
    fugacity_model: Literal[
        'vdW', 'PR', 'RK', 'SRK'
    ] = 'SRK',
    activity_model: Literal[
        'NRTL', 'UNIQUAC'
    ] = 'NRTL',
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
    solver_method: Literal[
        'root', 'least-squares', 'fsolve'
    ] = 'root',
    message: str | None = None,
    **kwargs
):
    '''
    The `bubble-point temperature` (BPT) calculation determines the temperature at which the first bubble of vapor forms when a liquid mixture is heated at a constant pressure. It helps identify the temperature at which the liquid will start vaporizing.

    Parameters
    ----------
    components : List[Component]
        List of Component objects representing the mixture.
    pressure : Pressure
        Pressure object specifying the pressure at which to calculate the bubble-point temperature.
    model_source : ModelSource
        ModelSource object containing the data source and equation of state information.
    equilibrium_model : Literal['raoult', 'modified-raoult'], optional
        The equilibrium model to use for the calculation. Options are 'raoult' or 'modified-raoult'. Default is 'raoult'.
    fugacity_model : Literal['vdW', 'PR', 'RK', 'SRK']]
        The fugacity model to use for the calculation. Options are 'vdW' (Van der Waals), 'PR' (Peng-Robinson), 'RK' (Redlich-Kwong), or 'SRK' (Soave-Redlich-Kwong). Default is SRK.
    activity_model : Literal['NRTL', 'UNIQUAC']
        The activity coefficient model to use for the calculation. Options are 'NRTL' (Non-Random Two-Liquid) or 'UNIQUAC' (Universal Quasi-Chemical). Default is NRTL.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State',
        'Formula-Name-State'], optional
        The key to identify components. Options are "Name-State" or "Formula-State". Default is "Name-State".
    solver_method : Literal['root', 'least-squares', 'fsolve'], optional
        The numerical solver method to use for the calculation. Options are 'root', 'least-squares', or 'fsolve'. Default is 'root'.
    message : str | None, optional
        Custom message to display during the calculation. Default is None.
    **kwargs : dict, optional
        Additional keyword arguments specific to the calculation.
        - activity_inputs : dict, model parameters.
        - activity_coefficients : dict, activity coefficients for each component.
        - guess_temperature : float, initial guess for temperature [K], default is 295 K.

    Returns
    -------
    res : dict
        Dictionary containing the results of the calculation.
        - bubble_temperature: bubble temperature [K]
        - pressure: pressure [Pa]
        - feed_mole_fraction: liquid mole fraction (xi)
        - vapor_mole_fraction: vapor mole fraction (yi)
        - liquid_mole_fraction: liquid mole fraction (xi)
        - vapor_pressure: vapor pressure [Pa]
        - activity_coefficient: activity coefficient [dimensionless]
        - K_ratio: K-ratio [dimensionless]
        - equilibrium_model: equilibrium model used for the calculation
        - fugacity_model: fugacity model used for the calculation
        - activity_model: activity coefficient model used for the calculation
        - components: list of components used in the calculation
        - message: message displayed during the calculation
        - solver_method: solver method used for the calculation
        - computation_time: computation time [s]

    Notes
    -----
    The summary of the calculation is as follows:

    - Known Information: Pressure (P) and mole fraction of the components in the liquid phase (xi).
    then zi = xi.
    - Computed Information: Temperature (T) and mole fraction in the vapor phase (yi).
    '''
    try:
        # SECTION: generate inputs
        input_dict = _input_generator(
            model_type='bpt',
            components=components,
            model_source=model_source,
            pressure=pressure,
            component_key=component_key,
        )

        # unpack inputs
        vle_model: VLE = input_dict['vle_model']
        inputs: Dict[str, Any] = input_dict['inputs']

        # SECTION: perform BPT calculation
        try:
            res = vle_model.bubble_temperature(
                inputs=inputs,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                solver_method=solver_method,
                message=message,
                **kwargs
            )
            return res
        except Exception as e:
            logger.error(f"BPT calculation failed!, {e}")
            raise

    except Exception as e:
        logger.error(f"Error in calc_bpt: {e}")
        raise e

# SECTION: Dew-point temperature (DPT) calculation


def calc_dew_point_temperature(
    components: List[Component],
    pressure: Pressure,
    model_source: ModelSource,
    equilibrium_model: Literal[
        'raoult',
        'modified-raoult'
    ] = 'raoult',
    fugacity_model: Literal[
        'vdW', 'PR', 'RK', 'SRK'
    ] = 'SRK',
    activity_model: Literal[
        'NRTL', 'UNIQUAC'
    ] = 'NRTL',
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
    solver_method: Literal[
        'root', 'least-squares', 'fsolve'
    ] = 'root',
    message: str | None = None,
    **kwargs
):
    '''
    The `dew-point temperature` (DPT) calculation determines the temperature at which the first drop of liquid condenses when a vapor mixture is cooled at a constant pressure. It identifies the temperature at which vapor will start to condense.

    Parameters
    ----------
    components : List[Component]
        List of Component objects representing the mixture.
    pressure : Pressure
        Pressure object specifying the pressure at which to calculate the bubble-point temperature.
    model_source : ModelSource
        ModelSource object containing the data source and equation of state information.
    equilibrium_model : Literal['raoult', 'modified-raoult'], optional
        The equilibrium model to use for the calculation. Options are 'raoult' or 'modified-raoult'. Default is 'raoult'.
    fugacity_model : Literal['vdW', 'PR', 'RK', 'SRK']]
        The fugacity model to use for the calculation. Options are 'vdW' (Van der Waals), 'PR' (Peng-Robinson), 'RK' (Redlich-Kwong), or 'SRK' (Soave-Redlich-Kwong). Default is SRK.
    activity_model : Literal['NRTL', 'UNIQUAC']
        The activity coefficient model to use for the calculation. Options are 'NRTL' (Non-Random Two-Liquid) or 'UNIQUAC' (Universal Quasi-Chemical). Default is NRTL.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State',
        'Formula-Name-State'], optional
        The key to identify components. Options are "Name-State" or "Formula-State". Default is "Name-State".
    solver_method : Literal['root', 'least-squares', 'fsolve'], optional
        The numerical solver method to use for the calculation. Options are 'root', 'least-squares', or 'fsolve'. Default is 'root'.
    message : str | None, optional
        Custom message to display during the calculation. Default is None.

    **kwargs : dict, optional
        Additional keyword arguments specific to the calculation.
        - activity_inputs : dict, model parameters.
        - activity_coefficients : dict, activity coefficients for each component.
        - guess_temperature : float, initial guess for temperature [K], default is 295 K.

    Returns
    -------
    res : dict
        Dictionary containing the results of the calculation.
        - dew_temperature: dew temperature [K]
        - pressure: pressure [Pa]
        - feed_mole_fraction: vapor mole fraction (yi)
        - liquid_mole_fraction: liquid mole fraction (xi)
        - vapor_pressure: vapor pressure [Pa]
        - activity_coefficient: activity coefficient [dimensionless]
        - K_ratio: K-ratio [dimensionless]
        - equilibrium_model: equilibrium model used for the calculation
        - fugacity_model: fugacity model used for the calculation
        - activity_model: activity coefficient model used for the calculation
        - components: list of components used in the calculation
        - message: message displayed during the calculation
        - solver_method: solver method used for the calculation
        - computation_time: computation time [s]

    Notes
    -----
    The summary of the calculation is as follows:

    - Known Information: Pressure (P) and mole fraction of the components in the vapor phase (yi).
    then zi = yi.
    - Computed Information: Temperature (T) and mole fraction in the liquid phase (xi).
    '''
    try:
        # SECTION: generate inputs
        input_dict = _input_generator(
            model_type='dpt',
            components=components,
            model_source=model_source,
            pressure=pressure,
            component_key=component_key,
        )

        # unpack inputs
        vle_model: VLE = input_dict['vle_model']
        inputs: Dict[str, Any] = input_dict['inputs']

        # SECTION: perform BPT calculation
        try:
            res = vle_model.dew_temperature(
                inputs=inputs,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                solver_method=solver_method,
                message=message,
                **kwargs
            )
            return res
        except Exception as e:
            logger.error(f"DPT calculation failed!, {e}")
            raise
    except Exception as e:
        logger.error(f"Error in calc_dpt: {e}")
        raise e

# SECTION: Isothermal Flash (IF) calculation


def calc_isothermal_flash(
    components: List[Component],
    temperature: Temperature,
    pressure: Pressure,
    model_source: ModelSource,
    equilibrium_model: Literal[
        'raoult', 'modified-raoult'
    ] = 'raoult',
    fugacity_model: Literal[
        'vdW', 'PR', 'RK', 'SRK'
    ] = 'SRK',
    activity_model: Literal[
        'NRTL', 'UNIQUAC'
    ] = 'NRTL',
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = "Name-State",
    solver_method: Literal[
        'least_squares', 'minimize'
    ] = 'least_squares',
    flash_checker: bool = False,
    message: str | None = None,
    **kwargs
) -> Dict[str, Any]:
    '''
    The `isothermal-flash` (IF) calculation This calculation determines the vapor and liquid phase compositions and amounts at a specified temperature and pressure. The system is "flashed" isothermally, meaning the temperature is kept constant while the phase behavior is calculated for a mixture.

    Parameters
    ----------
    components : List[Component]
        List of Component objects representing the mixture.
    temperature : Temperature
        Temperature object specifying the temperature for the flash calculation.
    pressure : Pressure
        Pressure object specifying the pressure for the flash calculation.
    model_source : ModelSource
        ModelSource object containing the data source and equation of state information.
    equilibrium_model : Literal['raoult', 'modified-raoult'], optional
        The equilibrium model to use for the calculation. Options are 'raoult' or 'modified-raoult'. Default is 'raoult'.
    fugacity_model : Literal['vdW', 'PR', 'RK', 'SRK']]
        The fugacity model to use for the calculation. Options are 'vdW' (Van der Waals), 'PR' (Peng-Robinson), 'RK' (Redlich-Kwong), or 'SRK' (Soave-Redlich-Kwong). Default is SRK.
    activity_model : Literal['NRTL', 'UNIQUAC']
        The activity coefficient model to use for the calculation. Options are 'NRTL' (Non-Random Two-Liquid) or 'UNIQUAC' (Universal Quasi-Chemical). Default is NRTL.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State',
        'Formula-Name-State'], optional
        The key to identify components. Options are "Name-State" or "Formula-State". Default is "Name-State".
    solver_method : Literal['least_squares', 'minimize'], optional
        The numerical solver method to use for the calculation. Options are 'least_squares' or
        'minimize'. Default is 'least_squares'.
    flash_checker : bool, optional
        If True, performs a flash calculation check before the main calculation. Default is False.
    message : str | None, optional
        Custom message to display during the calculation. Default is None.
    **kwargs : dict, optional
        Additional keyword arguments specific to the calculation.
        - activity_inputs : dict, model parameters.
        - activity_coefficients : dict, activity coefficients for each component.

    Returns
    -------
    res : dict
        Dictionary containing the results of the calculation.
        - vapor_to_liquid_ratio: V/F ratio [dimensionless]
        - liquid_to_vapor_ratio: L/F ratio [dimensionless]
        - feed_mole_fraction: liquid mole fraction (zi)
        - liquid_mole_fraction: liquid mole fraction (xi)
        - vapor_mole_fraction: vapor mole fraction (yi)
        - vapor_pressure: vapor pressure [Pa]
        - K_ratio: K-ratio [dimensionless]
        - temperature: temperature [K]
        - pressure: pressure [Pa]
        - equilibrium_model: equilibrium model used for the calculation
        - fugacity_model: fugacity model used for the calculation
        - activity_model: activity coefficient model used for the calculation
        - components: list of components used in the calculation
        - message: message displayed during the calculation
        - solver_method: solver method used for the calculation
        - flash_checker: flash checker used for the calculation
        - flash_checker_res: flash checker result
        - computation_time: computation time [s]

    Notes
    -----
    The summary of the calculation is as follows:

    - Known Information: Temperature (T), Pressure (P), and feed mole fraction of the components (zi).
    - Computed Information: Vapor and liquid phase compositions (yi and xi), and the vapor-to-liquid ratio (V/F).

    Flash occurs when the bubble pressure is greater than the dew pressure, and the system is in a two-phase region. The calculation will return the phase compositions and flow rates based on the specified temperature and pressure.

    flash case:
    - P[bubble] > P[flash] > P[dew] results in the two phase feed
    - P[bubble] < P[flash] results in the liquid phase feed
    - P[dew] > P[flash] results in the vapor phase feed
    '''
    try:
        # SECTION: generate inputs
        input_dict = _input_generator(
            model_type='if',
            components=components,
            model_source=model_source,
            temperature=temperature,
            pressure=pressure,
            component_key=component_key,
        )

        # unpack inputs
        vle_model: VLE = input_dict['vle_model']
        inputs: Dict[str, Any] = input_dict['inputs']

        # SECTION: perform IF calculation
        try:
            res = vle_model.flash_isothermal(
                inputs=inputs,
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model,
                solver_method=solver_method,
                flash_checker=flash_checker,
                message=message,
                **kwargs
            )
            return res
        except Exception as e:
            logger.error(f"IF calculation failed!, {e}")
            raise
    except Exception as e:
        logger.error(f"Error in calc_if: {e}")
        raise e

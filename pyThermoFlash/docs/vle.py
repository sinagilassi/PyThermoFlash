# VLE
# -----
# packages/modules
import time  # Import the time module
from typing import List, Dict, Optional, Literal, Any
import numpy as np
import pycuc
# local
from .equilibria import Equilibria
from .source import Source


class VLE(Equilibria):
    '''
    Model for vapor-liquid equilibrium (VLE) calculations for multi-component systems using the Raoult's law, modified Raoult's law and more.
    '''

    def __init__(
        self,
        components: List[str],
        model_source: Optional[Dict] = None,
        **kwargs
    ):
        '''
        Initialize the Partition class.
        '''
        self.__model_source = model_source
        self.__datasource = {
        } if model_source is None else model_source['datasource']
        self.__equationsource = {
        } if model_source is None else model_source['equationsource']

        # NOTE: init class
        Equilibria.__init__(self, components, model_source)

    @property
    def model_source(self) -> Dict:
        '''
        Get the model source property.

        Returns
        -------
        dict
            The model source dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__model_source is None:
            return {}
        return self.__model_source

    @property
    def datasource(self) -> Dict:
        '''
        Get the datasource property.

        Returns
        -------
        dict
            The datasource dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__datasource is None:
            return {}
        return self.__datasource

    @property
    def equationsource(self) -> Dict:
        '''
        Get the equationsource property.

        Returns
        -------
        dict
            The equationsource dictionary.
        '''
        # NOTE: check if model source is valid
        if self.__equationsource is None:
            return {}
        return self.__equationsource

    def __set_models(
        self,
        equilibrium_model: str,
        fugacity_model: Optional[str],
        activity_model: Optional[str]
    ):
        '''
        Set the models for the calculation.

        Parameters
        ----------
        equilibrium_model : str
            The equilibrium model to use for the calculation.
        fugacity_model : str, optional
            The fugacity model to use for the calculation.
        activity_model : str, optional
            The activity coefficient model to use for the calculation.

        Returns
        -------
        dict
            A dictionary containing the models used for the calculation.
            - equilibrium_model : str
                The equilibrium model used for the calculation.
            - fugacity_model : str, optional
                The fugacity model used for the calculation.
            - activity_model : str, optional
                The activity coefficient model used for the calculation.
        '''
        try:
            # NOTE: set based on equilibrium model
            if equilibrium_model == 'raoult':
                fugacity_model_ = None
                activity_model_ = None
            elif equilibrium_model == 'modified-raoult':
                fugacity_model_ = None
                activity_model_ = activity_model
            elif equilibrium_model == 'fugacity-ratio':
                fugacity_model_ = fugacity_model
                activity_model_ = None
            else:
                raise ValueError(
                    f"Invalid equilibrium model: {equilibrium_model}.")

            # res
            return {
                "equilibrium_model": equilibrium_model,
                "fugacity_model": fugacity_model_,
                "activity_model": activity_model_
            }

        except Exception as e:
            raise Exception(
                f"Error in setting models: {e}")

    def bubble_pressure(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        fugacity_model: Optional[Literal[
            'vdW', 'PR', 'RK', 'SRK'
        ]] = None,
        activity_model: Optional[Literal[
            'NRTL', 'UNIQUAC'
        ]] = None,
        message: Optional[str] = None,
        **kwargs
    ):
        '''
        The Bubble-Pressure (BP) calculation determines the pressure at which the first bubble of vapor forms when a liquid mixture is heated at a constant temperature. It is used to find the pressure for a given temperature at which the liquid will begin to vaporize.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Temperature at which to calculate the bubble pressure (in Kelvin, Celsius, Fahrenheit).
        equilibrium_model : str, optional
            The equilibrium model to use for the calculation. Default is 'raoult'.
        fugacity_model : str, optional
            The fugacity model to use for the calculation. Default is 'SRK'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
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


        For each component, the model parameters can be provided for activity models like NRTL or UNIQUAC:

        ```python
        # activity inputs
        activity_inputs = {
            'alpha': alpha,
            'a_ij': a_ij,
            'b_ij': b_ij,
            'c_ij': c_ij,
            'd_ij': d_ij
        }
        ```

        Or the activity coefficients can be provided directly:

        # calculated activity coefficients (for bubble-pressure calculation)
        ```python
        activity_coefficients = {
            'benzene': 1.0,
            'toluene': 1.1
        }
        ```
        '''
        try:
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'temperature' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction' and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: temperature [K]
            temperature = pycuc.convert_from_to(
                inputs['temperature'][0], inputs['temperature'][1], 'K')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check temperature
            if not isinstance(temperature, (int, float)):
                raise ValueError("Temperature must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                float(mole_fractions[component]) for component in components
            ]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = temperature

                # NOTE: execute
                _VaPr_res = VaPr_eq.cal(**_VaPr_args)
                # extract
                _VaPr_value = _VaPr_res['value']
                _VaPr_unit = _VaPr_res['unit']
                # unit conversion
                # NOTE: unit conversion
                _unit_block = f"{_VaPr_unit} => Pa"
                _VaPr = pycuc.to(_VaPr_value, _unit_block)
                # set
                VaPr_comp[component] = {
                    "value": _VaPr,
                    "unit": "Pa"}

            # SECTION: set models
            # check equilibrium model
            res_ = self.__set_models(
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model
            )

            # NOTE: parameters
            params = {
                "components": components,
                "mole_fraction": np.array(mole_fractions),
                "temperature": {
                    "value": temperature,
                    "unit": "K"
                },
                "vapor_pressure": VaPr_comp,
                "equilibrium_model": equilibrium_model,
                "fugacity_model": res_['fugacity_model'],
                "activity_model": res_['activity_model'],
            }

            # res
            res = self._BP(params, **kwargs)

            # NOTE: set message
            message = message if message is not None else "Bubble Pressure Calculation"
            # add
            res['message'] = message

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["fugacity_model"] = res_['fugacity_model']
            res["activity_model"] = res_['activity_model']

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # returns
            return res
        except Exception as e:
            raise ValueError(
                f"Error in bubble_pressure calculation: {e}")

    def dew_pressure(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        fugacity_model: Optional[Literal[
            'vdW', 'PR', 'RK', 'SRK'
        ]] = None,
        activity_model: Optional[Literal[
            'NRTL', 'UNIQUAC'
        ]] = None,
        message: Optional[str] = None,
        **kwargs
    ):
        '''
        The dew-point pressure (DP) calculation determines the pressure at which the first drop of liquid condenses when a vapor mixture is cooled at a constant temperature. It is used to find the pressure at which vapor will begin to condense.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Temperature at which to calculate the bubble pressure (in Kelvin, Celsius, Fahrenheit).
        equilibrium_model : str, optional
            The equilibrium model to use for the calculation. Default is 'raoult'.
        fugacity_model : str, optional
            The fugacity model to use for the calculation. Default is 'SRK'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
            - activity_inputs : dict, model parameters.
            - activity_coefficients : dict, activity coefficients for each component.

        Returns
        -------
        res : dict
            Dictionary containing the results of the calculation.
            - dew_pressure: dew pressure [Pa]
            - temperature: temperature [K]
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
            - computation_time: computation time [s]

        Notes
        -----
        - Temperature must be in Kelvin, Celsius, or Fahrenheit, and will be converted to Kelvin.
        - Mole fractions must be between 0 and 1.
        - The function will raise a ValueError if the inputs are not valid.
        - The default equilibrium model is 'raoult', and the default activity model is 'NRTL'.
        - The function will return the dew pressure in Pascals.
        '''
        try:
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'temperature' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction' and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: temperature [K]
            temperature = pycuc.convert_from_to(
                inputs['temperature'][0], inputs['temperature'][1], 'K')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check temperature
            if not isinstance(temperature, (int, float)):
                raise ValueError("Temperature must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                mole_fractions[component] for component in components
            ]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = temperature

                # NOTE: execute
                _VaPr_res = VaPr_eq.cal(**_VaPr_args)
                # extract
                _VaPr_value = _VaPr_res['value']
                _VaPr_unit = _VaPr_res['unit']
                # unit conversion
                # NOTE: unit conversion
                _unit_block = f"{_VaPr_unit} => Pa"
                _VaPr = pycuc.to(_VaPr_value, _unit_block)
                # set
                VaPr_comp[component] = {
                    "value": _VaPr,
                    "unit": "Pa",
                    "equation": VaPr_eq,
                    "args": _VaPr_args,
                    "return": VaPr_eq.returns
                }

            # SECTION: set models
            # check equilibrium model
            res_ = self.__set_models(
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model
            )

            # NOTE: parameters
            params = {
                "components": components,
                "mole_fraction": np.array(mole_fractions),
                "temperature": {
                    "value": temperature,
                    "unit": "K"
                },
                "vapor_pressure": VaPr_comp,
                "equilibrium_model": equilibrium_model,
                "fugacity_model": res_['fugacity_model'],
                "activity_model": res_['activity_model'],
            }

            # res
            res = self._DP(params, **kwargs)

            # NOTE: set message
            message = message if message is not None else "Dew Pressure Calculation"
            # add
            res['message'] = message

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["fugacity_model"] = res_['fugacity_model']
            res["activity_model"] = res_['activity_model']

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # returns
            return res
        except Exception as e:
            raise ValueError(
                f"Error in bubble_pressure calculation: {e}")

    def bubble_temperature(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        fugacity_model: Literal[
            'vdW', 'PR', 'RK', 'SRK'
        ] = 'SRK',
        activity_model: Literal[
            'NRTL', 'UNIQUAC']
        = 'NRTL',
        solver_method: Literal[
            'root', 'least-squares', 'fsolve'
        ] = 'root',
        message: Optional[str] = None,
        **kwargs
    ):
        '''
        The `bubble-point temperature` (BT) calculation determines the temperature at which the first bubble of vapor forms when a liquid mixture is heated at a constant pressure. It helps identify the temperature at which the liquid will start vaporizing.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Temperature at which to calculate the bubble pressure (in Kelvin, Celsius, Fahrenheit).
        equilibrium_model : str, optional
            The equilibrium model to use for the calculation. Default is 'raoult'.
        fugacity_model : str, optional
            The fugacity model to use for the calculation. Default is 'SRK'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        solver_method : str, optional
            The solver method to use for the calculation. Default is 'root'.
            - 'root' : Root-finding algorithm.
            - 'least-squares' : Least-squares optimization algorithm.
            - 'fsolve' : SciPy's fsolve function.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
            - `guess_temperature` : float
                Initial guess for the bubble-point temperature (in Kelvin, Celsius, Fahrenheit), default is 295 K.
            - activity_inputs : dict,
                model parameters.
            - activity_coefficients : dict,
                activity coefficients for each component.

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
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'pressure' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction' and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: pressure [Pa]
            pressure = pycuc.convert_from_to(
                inputs['pressure'][0], inputs['pressure'][1], 'Pa')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check temperature
            if not isinstance(pressure, (int, float)):
                raise ValueError("pressure must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                mole_fractions[component] for component in components
            ]

            # SECTION: check for activity model

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                VaPr_eq = None
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = None

                # set
                VaPr_comp[component] = {
                    "value": VaPr_eq,
                    "args": _VaPr_args,
                    "return": VaPr_eq.returns
                }

            # SECTION: set models
            # check equilibrium model
            res_ = self.__set_models(
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model
            )

            # NOTE: parameters
            params = {
                "components": components,
                "mole_fraction": np.array(mole_fractions),
                "pressure": {
                    "value": pressure,
                    "unit": "Pa"
                },
                "vapor_pressure": VaPr_comp,
                "equilibrium_model": equilibrium_model,
                "fugacity_model": res_['fugacity_model'],
                "activity_model": res_['activity_model'],
                "solver_method": solver_method
            }

            # res
            res = self._BT(params, **kwargs)

            # NOTE: set message
            message = message if message is not None else "Bubble Temperature Calculation"
            # add
            res['message'] = message

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["fugacity_model"] = res_['fugacity_model']
            res["activity_model"] = res_['activity_model']
            res["solver_method"] = solver_method

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # returns
            return res
        except Exception as e:
            raise Exception(
                f"Error in bubble_temperature calculation: {e}")

    def dew_temperature(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        fugacity_model: Literal[
            'vdW', 'PR', 'RK', 'SRK'
        ] = 'SRK',
        activity_model: Literal[
            'NRTL', 'UNIQUAC']
        = 'NRTL',
        solver_method: Literal[
            'root', 'least-squares', 'fsolve'
        ] = 'least-squares',
        message: Optional[str] = None,
            **kwargs
    ):
        '''
        The `dew-point temperature` (DT) calculation determines the temperature at which the first drop of liquid condenses when a vapor mixture is cooled at a constant pressure. It identifies the temperature at which vapor will start to condense.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Temperature at which to calculate the bubble pressure (in Kelvin, Celsius, Fahrenheit).
        equilibrium_model : str, optional
            The equilibrium model to use for the calculation. Default is 'raoult'.
        fugacity_model : str, optional
            The fugacity model to use for the calculation. Default is 'SRK'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        solver_method : str, optional
            The solver method to use for the calculation. Default is 'root'.
            - 'root' : Root-finding algorithm.
            - 'least-squares' : Least-squares optimization algorithm.
            - 'fsolve' : SciPy's fsolve function.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
            - `guess_temperature` : float
                Initial guess for the bubble-point temperature (in Kelvin, Celsius, Fahrenheit), default is 295 K.

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
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'pressure' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction' and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: pressure [Pa]
            # set
            pressure_value = float(inputs['pressure'][0])
            pressure_unit = str(inputs['pressure'][1])

            pressure = pycuc.convert_from_to(
                pressure_value, pressure_unit, 'Pa')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check temperature
            if not isinstance(pressure, (int, float)):
                raise ValueError("pressure must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                mole_fractions[component] for component in components
            ]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                VaPr_eq = None
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = None

                # set
                VaPr_comp[component] = {
                    "value": VaPr_eq,
                    "args": _VaPr_args,
                    "return": VaPr_eq.returns
                }

            # SECTION: set models
            # check equilibrium model
            res_ = self.__set_models(
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model
            )

            # NOTE: parameters
            params = {
                "components": components,
                "mole_fraction": np.array(mole_fractions),
                "pressure": {
                    "value": pressure,
                    "unit": "Pa"
                },
                "vapor_pressure": VaPr_comp,
                "equilibrium_model": equilibrium_model,
                "fugacity_model": res_['fugacity_model'],
                "activity_model": res_['activity_model'],
                "solver_method": solver_method
            }

            # res
            res = self._DT(params, **kwargs)

            # NOTE: set message
            message = message if message is not None else "Dew Temperature Calculation"
            # add
            res['message'] = message

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["fugacity_model"] = res_['fugacity_model']
            res["activity_model"] = res_['activity_model']
            res["solver_method"] = solver_method

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # returns
            return res
        except Exception as e:
            raise Exception(
                f"Error in dew_temperature calculation: {e}")

    def flash_isothermal(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        fugacity_model: Literal[
            'vdW', 'PR', 'RK', 'SRK'
        ] = 'SRK',
        activity_model: Literal[
            'NRTL', 'UNIQUAC']
        = 'NRTL',
        solver_method: Literal[
            'least_squares', 'minimize'
        ] = 'least_squares',
        flash_checker: bool = False,
        message: Optional[str] = None,
        **kwargs
    ):
        '''
        The `isothermal-flash` (IFL) calculation This calculation determines the vapor and liquid phase compositions and amounts at a specified temperature and pressure. The system is "flashed" isothermally, meaning the temperature is kept constant while the phase behavior is calculated for a mixture.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Flash temperature (in Kelvin, Celsius, Fahrenheit).
            - pressure : float
                Flash pressure (in Pascal, bar, atm).
        equilibrium_model : str, optional
            the equilibrium model to use for the calculation. Default is 'raoult'.
        fugacity_model : str, optional
            The fugacity model to use for the calculation. Default is 'SRK'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        solver_method : str, optional
            The solver method to use for the calculation. Default is 'least_squares'.
            - 'minimize' : Minimize the objective function.
            - 'least-squares' : Least-squares optimization algorithm.
        flash_checker : bool, optional
            If True, check if the system is in a two-phase region. Default is False.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
            - `guess_V_F_ratio`: initial guess for the vapor-to-liquid ratio (V/F), default is 0.5

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
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'pressure' not in inputs or 'temperature' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction', 'pressure', and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: pressure [Pa]
            pressure = pycuc.convert_from_to(
                inputs['pressure'][0], inputs['pressure'][1], 'Pa')
            # NOTE: temperature [K]
            temperature = pycuc.convert_from_to(
                inputs['temperature'][0], inputs['temperature'][1], 'K')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check pressure
            if not isinstance(pressure, (int, float)):
                raise ValueError("pressure must be numeric.")
            # check temperature
            if not isinstance(temperature, (int, float)):
                raise ValueError("temperature must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                mole_fractions[component] for component in components
            ]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = temperature

                # set
                VaPr_comp[component] = {
                    "value": VaPr_eq,
                    "args": _VaPr_args,
                    "return": VaPr_eq.returns
                }

            # SECTION: set models
            # check equilibrium model
            res_ = self.__set_models(
                equilibrium_model=equilibrium_model,
                fugacity_model=fugacity_model,
                activity_model=activity_model
            )

            # NOTE: parameters
            params = {
                "components": components,
                "mole_fraction": np.array(mole_fractions),
                "pressure": {
                    "value": pressure,
                    "unit": "Pa"
                },
                "temperature": {
                    "value": temperature,
                    "unit": "K"
                },
                "vapor_pressure": VaPr_comp,
                "equilibrium_model": equilibrium_model,
                "fugacity_model": res_['fugacity_model'],
                "activity_model": res_['activity_model'],
                "solver_method": solver_method
            }

            # SECTION: flash checker
            flash_checker_res_ = None
            if flash_checker:
                # check
                flash_checker_res_ = self._flash_checker(
                    z_i=mole_fractions,
                    Pf=pressure,
                    Tf=temperature,
                    VaPr_comp=VaPr_comp,
                )

                # check
                if flash_checker_res_ is False:
                    raise ValueError(
                        "Flash calculation failed! The system is not in a two-phase region.")

            # SECTION: flash calculation
            res = self._IFL(params, **kwargs)

            # NOTE: set message
            message = message if message is not None else "Flash Isothermal Calculation"
            # add
            res['message'] = message

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["flash_checker"] = flash_checker
            res["flash_checker_res"] = flash_checker_res_
            res["fugacity_model"] = res_['fugacity_model']
            res["activity_model"] = res_['activity_model']
            res["solver_method"] = solver_method

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # res
            return res
        except Exception as e:
            raise Exception(
                f"Error in flash_isothermal calculation: {e}")

    def is_flashable(
        self,
        inputs: Dict[str, Any],
        equilibrium_model: Literal[
            'raoult', 'modified-raoult'
        ] = 'raoult',
        message: Optional[str] = None,
        **kwargs
    ):
        '''
        Check if the mixture is flashable at the given pressure and temperature based on Raoult's law.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Flash temperature (in Kelvin, Celsius, Fahrenheit).
            - pressure : float
                Flash pressure (in Pascal, bar, atm).
        equilibrium_model : str, optional
            the equilibrium model to use for the calculation. Default is 'raoult'.
        message : str, optional
            Message to display during the calculation. Default is None.
        **kwargs : dict, optional
            Additional parameters for the model.
            - `guess_V_F_ratio`: initial guess for the vapor-to-liquid ratio (V/F), default is 0.5

        Returns
        -------
        res : dict
            Dictionary containing the results of the calculation.
            - feed_mole_fraction: feed mole fraction (zi)
            - vapor_pressure: vapor pressure [Pa]
            - temperature: temperature [K]
            - pressure: pressure [Pa]
            - equilibrium_model: equilibrium model used for the calculation
            - components: list of components used in the calculation
            - message: message displayed during the calculation
            - flash_checker_res: flash checker result
            - computation_time: computation time [s]

        Notes
        -----
        Flash occurs when the bubble pressure is greater than the dew pressure, and the system is in a two-phase region. The calculation will return the phase compositions and flow rates based on the specified temperature and pressure.

        flash case:
        - P[bubble] > P[flash] > P[dew] results in the two phase feed
        - P[bubble] < P[flash] results in the liquid phase feed
        - P[dew] > P[flash] results in the vapor phase feed
        '''
        try:
            # ! Start timing
            start_time = time.time()

            # SECTION: check inputs
            # check inputs
            if not isinstance(inputs, dict):
                raise ValueError("Inputs must be a dictionary.")

            # check keys
            if 'mole_fraction' not in inputs or 'pressure' not in inputs or 'temperature' not in inputs:
                raise ValueError(
                    "Inputs must contain 'mole_fraction', 'pressure', and 'temperature' keys.")

            # NOTE: set inputs - dict
            mole_fractions = inputs['mole_fraction']
            # NOTE: pressure [Pa]
            pressure = pycuc.convert_from_to(
                inputs['pressure'][0], inputs['pressure'][1], 'Pa')
            # NOTE: temperature [K]
            temperature = pycuc.convert_from_to(
                inputs['temperature'][0], inputs['temperature'][1], 'K')

            # check mole fractions
            if not isinstance(mole_fractions, dict):
                raise ValueError("Mole fractions must be a dictionary.")

            if not all(isinstance(v, (int, float)) for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be numeric.")

            if not all(0 <= v <= 1 for v in mole_fractions.values()):
                raise ValueError("Mole fractions must be between 0 and 1.")

            # check pressure
            if not isinstance(pressure, (int, float)):
                raise ValueError("pressure must be numeric.")
            # check temperature
            if not isinstance(temperature, (int, float)):
                raise ValueError("temperature must be numeric.")

            # NOTE: components
            components = self.components

            # mole fractions based on components id
            mole_fractions = [
                mole_fractions[component] for component in components
            ]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = Source_.eq_extractor(component, 'VaPr')

                # NOTE: args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(component, VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(component, VaPr_args_required)

                # NOTE: update P and T
                _VaPr_args['T'] = temperature

                # set
                VaPr_comp[component] = {
                    "value": VaPr_eq,
                    "args": _VaPr_args,
                    "return": VaPr_eq.returns
                }

            # SECTION: flash checker
            # check
            flash_checker_res = self._flash_checker(
                z_i=mole_fractions,
                Pf=pressure,
                Tf=temperature,
                VaPr_comp=VaPr_comp,
            )

            # NOTE: init res
            res = {}

            # NOTE: set message
            message = message if message is not None else "Is Flashable Calculation"
            # add
            res['message'] = message

            # NOTE: pressure and temperature
            res['pressure'] = {
                "value": pressure,
                "unit": "Pa"
            }
            res['temperature'] = {
                "value": temperature,
                "unit": "K"
            }

            # NOTE: mole fractions
            res['feed_mole_fraction'] = mole_fractions

            # NOTE: components
            res['components'] = components

            # NOTE: models
            res["equilibrium_model"] = equilibrium_model
            res["flash_checker_res"] = flash_checker_res

            # NOTE: set time
            # ! Stop timing
            end_time = time.time()
            computation_time = end_time - start_time
            # add to res
            res['computation_time'] = {
                "value": computation_time,
                "unit": "s"
            }

            # res
            return res
        except Exception as e:
            raise Exception(
                f"Error in is_flashable calculation: {e}")

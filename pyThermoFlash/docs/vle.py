# POOL
# -----
# packages/modules
from typing import List, Dict, Optional, Literal
import numpy as np
import pycuc
# local
from .equilibria import Equilibria
from .source import Source


class VLE(Equilibria):
    '''
    Model for vapor-liquid equilibrium (VLE) calculations for multi-component systems using the Raoult's law, modified Raoult's law and more.
    '''

    def __init__(self,
                 components: List[str],
                 model_source: Optional[Dict] = None,
                 **kwargs):
        '''
        Initialize the Partition class.
        '''
        self.__model_source = model_source

        # init class
        Equilibria.__init__(self, components)

    def bubble_pressure(self,
                        inputs: Dict[str, float],
                        equilibrium_model: Literal[
                            'raoult', 'modified-raoult', 'fugacity-ratio'
                        ] = 'raoult',
                        activity_model: Literal[
                            'NRTL', 'UNIQUAC']
                        = 'NRTL',
                        **kwargs):
        '''
        The Bubble-Pressure (BP) calculation determines the pressure at which the first bubble of vapor forms when a liquid mixture is heated at a constant temperature. It is used to find the pressure for a given temperature at which the liquid will begin to vaporize.

        Parameters
        ----------
        inputs : dict
            Dictionary containing the input parameters for the calculation.
            - mole_fraction : dict
                Dictionary of component names and their respective mole fractions.
            - temperature : float
                Temperature at which to calculate the bubble pressure (in Kelvin).
        equilibrium_model : str, optional
            The equilibrium model to use for the calculation. Default is 'raoult'.
        activity_model : str, optional
            The activity coefficient model to use for the calculation. Default is 'NRTL'.
        **kwargs : dict, optional
            Additional parameters for the model.


        '''
        try:
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
            mole_fractions = [mole_fractions[component]
                              for component in components]

            # SECTION: extract source
            Source_ = Source(self.model_source)

            # NOTE: vapor pressure source
            VaPr_eq = {}
            VaPr_comp = {}

            # looping through components
            for component in components:
                # NOTE: equation source
                # antoine equations [Pa]
                VaPr_eq = Source_.data_extractor(component, 'VaPr')
                # args
                VaPr_args = VaPr_eq.args
                # check args (SI)
                VaPr_args_required = Source_.check_args(VaPr_args)

                # build args
                _VaPr_args = Source_.build_args(VaPr_args_required)

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
                _vapor_pressure.append(_VaPr)

            # NOTE: parameters
            params = {
                "mole_fraction": np.array(mole_fractions),
                "temperature": temperature,
                "vapor_pressure": {
                    "equation": VaPr_eq,
                },
                "equilibrium_model": equilibrium_model,
                "activity_model": activity_model
            }

            # res
            return self.bubblePressure(params)
        except Exception as e:
            raise ValueError(
                f"Error in bubble_pressure calculation: {e}")

    def dew_pressure(self):
        pass

    def bubble_temperature(self, mole_fractions, pressure, guess_temperature=350, vapor_pressure_method='polynomial'):
        '''
        bubble temperature calculation

        args:
            mole_fraction: feed mole fraction (zi=xi)
            pressure: system pressure [Pa]
        '''
        # params
        params = {
            "zi": np.array(mole_fractions),
            "P": pressure
        }

        # config
        config = {
            "Tg0": guess_temperature,
            "VaPeCal": vapor_pressure_method
        }

        # cal
        _res = self.bubbleTemperature(params, config)

        # res
        return _res

    def dew_temperature(self, mole_fractions, pressure, guess_temperature=350, vapor_pressure_method='polynomial'):
        '''
        bubble temperature calculation

        args:
            mole_fraction: feed mole fraction (zi=xi)
            pressure: system pressure [Pa]
        '''
        # params
        params = {
            "zi": np.array(mole_fractions),
            "P": pressure
        }

        # config
        config = {
            "Tg0": guess_temperature,
            "VaPeCal": vapor_pressure_method
        }

        # cal
        _res = self.dewTemperature(params, config)

        # res
        return _res

    def flash_isothermal(self, mole_fractions, flash_pressure, flash_temperature, feed_pressure, feed_flowrate=1, guess_V_F_ratio=0.5, vapor_pressure_method='polynomial', model="raoult", activity_coefficient_model='van-laar'):
        '''
        isothermal flash calculation

        args:
            mole_fractions: feed mole fraction [-]
            flash_pressure: flash pressure [Pa]
            flash_temperature: flash temperature [K] - isothermal condition T[in]=T[out]
            feed_pressure: feed pressure to check the current state of the feed
            feed_flowrate: feed flowrate (mole basis) [mol/s]
            guess_V_F_ratio:
            vapor_pressure_method:
            model: "raoult"

        notes:
            flash case: P[bubble]>P[flash]>P[dew]
            cases:
                1. P[bubble]<P[flash] results in the liquid phase feed
                2. P[dew]>P[flash] results in the vapor phase feed
        '''
        # vapor pressure at inlet temperature
        VaPr = self.vaporPressureMixture(
            flash_temperature, vapor_pressure_method)

        # bubble pressure [Pa]
        BuPr = self.calBubblePressure(mole_fractions, VaPr)

        # dew pressure
        DePr = self.calDewPressure(mole_fractions, VaPr)

        # check flash
        flashState = False
        if BuPr > flash_pressure and flash_pressure > DePr:
            # two phase exists
            flashState = True
        else:
            # one phase exist (liquid)
            flashState = False

        # check current state of the feed
        feedState = True
        if feed_pressure < BuPr:
            feedState = False

        # params
        params = {
            "F": feed_flowrate,
            "zi": np.array(mole_fractions),
            "P_flash": flash_pressure,
            "T_flash": flash_temperature,
            "VaPe": VaPr
        }

        # config
        config = {
            "guess_V_F_ratio": guess_V_F_ratio,
            "VaPeCal": vapor_pressure_method,
            "model": model
        }

        # flash calculation
        V_F_ratio, L_F_ratio, xi, yi = self.flashIsothermalV2(params, config)

        # NOTE
        # ! display results
        # vapor pressure
        _display = []

        # header
        _display.append(['Parameter', 'Value', 'Unit'])

        for i in range(self.comp_num):
            _display.append(
                [self.components[i].symbol+' vapor-pressure [P*]', roundNum(VaPr[i], 3), 'Pa'])

        for i in range(self.comp_num):
            _display.append([self.components[i].symbol +
                            ' liquid mole fraction [x]', roundNum(xi[i], 4), '-'])

        for i in range(self.comp_num):
            _display.append([self.components[i].symbol +
                            ' vapor mole fraction [y]', roundNum(yi[i], 4), '-'])

        #
        _display.extend(
            [['bubble pressure', roundNum(BuPr, 2), 'Pa'],
             ['dew pressure', roundNum(DePr, 2), 'Pa'],
             ['feed flowrate', 1, 'mol/s'],
             ['liquid flowrate', roundNum(L_F_ratio, 3), 'mol/s'],
             ['vapor flowrate', roundNum(V_F_ratio, 3), 'mol/s']]
        )

        # log
        # self.colDisplay(_display)

        # res
        return flashState, BuPr, DePr, VaPr, V_F_ratio, L_F_ratio, xi, yi

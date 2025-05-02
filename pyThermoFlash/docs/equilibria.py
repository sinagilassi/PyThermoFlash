# VAPOR-LIQUID EQUILIBRIUM (VLE) CALCULATIONS
# --------------------------------------------
# import libs
from typing import List, Dict, Optional
import numpy as np
from scipy import optimize
import pycuc
import pyThermoModels as ptm
# local
from .source import Source
from .activity import Activity


class Equilibria:
    '''
    Phase Equilibria Calculations
    '''

    def __init__(self,
                 components: List[str],
                 **kwargs):
        '''Initialize the Equilibria class.'''
        components_ = [component.strip() for component in components]
        self.__components = components_
        self.__comp_num = len(components_)

    def __repr__(self):
        des = "Phase Equilibria Calculations"
        return des

    @property
    def components(self):
        '''Get the components of the system.'''
        return self.__components

    @property
    def component_num(self):
        '''Get the number of components in the system.'''
        return self.__comp_num

    def __mole_fraction_comp(self, z_i: List[float]) -> Dict[str, float]:
        '''
        Convert mole fraction list to dictionary.

        Parameters
        ----------
        z_i : list
            List of mole fractions.

        Returns
        -------
        dict
            Dictionary of mole fractions.
        '''
        # convert to dict
        return {str(self.components[i]): z_i[i]
                for i in range(self.component_num)}

    def _BP(self, params, **kwargs):
        '''
        The bubble-point pressure (BP) calculation determines the pressure at which the first bubble of vapor forms when a liquid mixture is heated at a constant temperature. It is used to find the pressure for a given temperature at which the liquid will begin to vaporize. This calculation is based on `Raoult's law`, which states that the vapor pressure of each component in the mixture is proportional to its mole fraction in the liquid phase.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : liquid mole fraction
            - T: temperature [K]
        kwargs : dict
            additional parameters for the calculation
                - `activity_inputs`: information for the calculation


        Returns
        -------
        dict
            Dictionary containing the following:
            - bubble_pressure: bubble pressure [Pa]
            - temperature: temperature [K]
            - feed_mole_fraction: liquid mole fraction (xi)
            - vapor_mole_fraction: vapor mole fraction (yi)
            - liquid_mole_fraction: liquid mole fraction (xi)
            - vapor_pressure: vapor pressure [Pa]
            - activity_coefficient: activity coefficient [dimensionless]
            - K_ratio: K-ratio [dimensionless]

        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Temperature (T) and mole fraction of the components in the liquid phase (xi).
        then zi = xi.
        - Computed Information: Pressure (P) and mole fraction in the vapor phase (yi).

        - The solution is obtained by the following steps:
            1. Calculate the vapor pressure of each component at the given temperature (T).
            2. Calculate the bubble pressure (P) using `Raoult's law` or other methods.
            3. Calculate the mole fraction in the vapor phase (yi) using `Raoult's law` or other methods.

        - inputs hold the information for the calculation:
            `alpha_ij`: binary interaction parameter
            `dg_ij`: dg_ij parameter
            `dU_ij`: dU_ij parameter
            `a_ij`: a_ij parameter
            `b_ij`: b_ij parameter
            `c_ij`: c_ij parameter
            `d_ij`: d_ij parameter
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            z_i = params['mole_fraction']
            # temperature [K]
            T = params['temperature']
            # equilibrium model
            eq_model = params['equilibrium_model']
            # fugacity model
            fugacity_model = params['fugacity_model']
            # activity model
            activity_model = params['activity_model']

            # NOTE: vapor pressure [Pa]
            VaPr_comp = params['vapor_pressure']

            # NOTE: mole fraction dict
            # convert to dict
            z_i_comp = self.__mole_fraction_comp(z_i)

            # temperature [K]
            T_value = T['value']

            # SECTION: activity coefficient
            # NOTE: init model
            # init NRTL model
            activity = Activity()

            # NOTE: check model
            if activity_model == 'NRTL':
                # calculate activity
                res_ = activity.NRTL(
                    self.components, z_i_comp, T_value, **kwargs)
                # extract
                AcCo_i = res_['value']
            elif activity_model == 'UNIQUAC':
                # calculate activity
                res_ = activity.UNIQUAC(
                    self.components, z_i_comp, T_value, **kwargs)
                # extract
                AcCo_i = res_['value']
            else:
                # equals unity for ideal solution
                AcCo_i = np.ones(self.component_num)

            # SECTION: vapor pressure calculation
            # vapor pressure [Pa]
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [Pa]
                VaPr_i[i] = VaPr_comp[component]['value']

            # bubble pressure [Pa]
            BuPr = np.sum(AcCo_i*z_i*VaPr_i)

            # vapor mole fraction
            y_i = np.zeros(self.component_num)

            # looping over components
            for i in range(self.component_num):
                y_i[i] = z_i[i]*VaPr_i[i]*AcCo_i[i]/BuPr

            # NOTE: Ki ratio
            K_i = np.multiply(y_i, 1/z_i)

            # NOTE: results
            res = {
                "bubble_pressure": {
                    "value": float(BuPr),
                    "unit": "Pa"
                },
                "temperature": {
                    "value": T_value,
                    "unit": "K"
                },
                "feed_mole_fraction": z_i,
                "vapor_mole_fraction": y_i,
                "liquid_mole_fraction": z_i,
                "vapor_pressure": {
                    "value": VaPr_i,
                    "unit": "Pa"
                },
                "activity_coefficient": {
                    "value": AcCo_i,
                    "unit": "dimensionless"
                },
                "K_ratio": {
                    "value": K_i,
                    "unit": "dimensionless"
                }
            }

            # res
            return res
        except Exception as e:
            raise Exception(f'bubble pressure calculation failed! {e}')

    def _DP(self, params, **kwargs):
        '''
        The dew-point pressure (DP) calculation determines the pressure at which the first drop of liquid condenses when a vapor mixture is cooled at a constant temperature. It is used to find the pressure at which vapor will begin to condense.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : vapor mole fraction (yi)
            - T: temperature [K]
        kwargs : dict
            additional parameters for the calculation

        Returns
        -------
        res : dict
            Dictionary containing the following:
            - dew_pressure: dew pressure [Pa]
            - temperature: temperature [K]
            - feed_mole_fraction: vapor mole fraction (yi)
            - vapor_mole_fraction: vapor mole fraction (yi)
            - liquid_mole_fraction: liquid mole fraction (xi)
            - vapor_pressure: vapor pressure [Pa]
            - activity_coefficient: activity coefficient [dimensionless]
            - K_ratio: K-ratio [dimensionless]



        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Temperature (T) and mole fraction of the components in the vapor phase (yi).
        then zi = yi.
        - Computed Information: Pressure (P) and mole fraction in the liquid phase (xi).


        - The solution is obtained by the following steps:
            1. Calculate the vapor pressure of each component at the given temperature (T).
            2. Calculate the dew pressure (P) using Raoult's law.
            3. Calculate the mole fraction in the liquid phase (xi) using `Raoult's law`.
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            y_i = params['mole_fraction']
            # temperature [K]
            T = params['temperature']
            # equilibrium model
            eq_model = params['equilibrium_model']
            # fugacity model
            fugacity_model = params['fugacity_model']
            # activity model
            activity_model = params['activity_model']

            # SECTION: activity coefficient
            # NOTE: mole fraction dict
            # convert to dict
            y_i_comp = self.__mole_fraction_comp(y_i)

            # temperature [K]
            T_value = T['value']

            # NOTE: init model
            # init NRTL model
            activity = Activity()

            # NOTE: check model
            if activity_model == 'NRTL':
                # calculate activity
                res_ = activity.NRTL(
                    self.components, y_i_comp, T_value, **kwargs)
                # extract
                AcCo_i = res_['value']
            elif activity_model == 'UNIQUAC':
                # calculate activity
                res_ = activity.UNIQUAC(
                    self.components, y_i_comp, T_value, **kwargs)
                # extract
                AcCo_i = res_['value']
            else:
                # equals unity for ideal solution
                AcCo_i = np.ones(self.component_num)

            # NOTE: vapor pressure [Pa]
            VaPr_comp = params['vapor_pressure']

            # vapor pressure [Pa]
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [Pa]
                VaPr_i[i] = VaPr_comp[component]['value']

            # NOTE: dew pressure [Pa]
            DePr = 1/np.dot(y_i, 1/(VaPr_i*AcCo_i))

            # NOTE: liquid mole fraction
            x_i = np.zeros(self.component_num)

            # looping over components
            for i in range(self.component_num):
                # mole fraction
                x_i[i] = y_i[i]*DePr/VaPr_i[i]*AcCo_i[i]

            # NOTE: k-ratio
            K_i = np.multiply(y_i, 1/x_i)

            # res
            res = {
                "dew_pressure": {
                    "value": float(DePr),
                    "unit": "Pa"
                },
                "temperature": {
                    "value": T,
                    "unit": "K"
                },
                "feed_mole_fraction": y_i,
                "vapor_mole_fraction": y_i,
                "liquid_mole_fraction": x_i,
                "vapor_pressure": {
                    "value": VaPr_i,
                    "unit": "Pa"
                },
                "activity_coefficient": {
                    "value": AcCo_i,
                    "unit": "dimensionless"
                },
                "K_ratio": {
                    "value": K_i,
                    "unit": "dimensionless"
                }
            }

            # res
            return res
        except Exception as e:
            raise Exception(f'dew pressure calculation failed! {e}')

    def _BT(self, params, **kwargs):
        '''
        The `bubble-point temperature` (BT) calculation determines the temperature at which the first bubble of vapor forms when a liquid mixture is heated at a constant pressure. It helps identify the temperature at which the liquid will start vaporizing.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : liquid mole fraction (xi)
            - P: pressure [Pa]
        kwargs : dict
            additional parameters for the calculation
            - `guess_temperature`: initial guess temperature [K], default is 295 K
            - `activity_inputs`: activity model inputs consists of:
                - activity_model: activity model name (NRTL, UNIQUAC)
                - alpha_ij: binary interaction parameter
                - dg_ij: dg_ij parameter
                - dU_ij: dU_ij parameter
                - a_ij: a_ij parameter
                - b_ij: b_ij parameter
                - c_ij: c_ij parameter
                - d_ij: d_ij parameter

        Returns
        -------
        res : dict
            Dictionary containing the following:
            - bubble_temperature: bubble temperature [K]
            - pressure: pressure [Pa]
            - feed_mole_fraction: liquid mole fraction (xi)
            - vapor_mole_fraction: vapor mole fraction (yi)
            - liquid_mole_fraction: liquid mole fraction (xi)
            - vapor_pressure: vapor pressure [Pa]
            - activity_coefficient: activity coefficient [dimensionless]
            - K_ratio: K-ratio [dimensionless]

        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Pressure (P) and mole fraction of the components in the liquid phase (xi).
        then zi = xi.
        - Computed Information: Temperature (T) and mole fraction in the vapor phase (yi).

        The solution is obtained by the following steps:
            1. Choose an initial guess for the bubble-point temperature (Tg).
            2. Calculate the vapor pressure of each component at the guessed temperature (Tg).
            3. Calculate the bubble pressure (Pb) using `Raoult's law` or other methods.
            4. Check if the calculated bubble pressure (Pb) equals the given pressure (P).
            5. Calculate the vapor mole fraction (yi) using `Raoult's law` or other methods.
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            z_i = params['mole_fraction']
            # pressure [Pa]
            P = params['pressure']
            # equilibrium model
            eq_model = params['equilibrium_model']
            # fugacity model
            fugacity_model = params['fugacity_model']
            # activity model
            activity_model = params['activity_model']
            # solver method
            solver_method = params['solver_method']

            # SECTION: vapor pressure calculation
            # NOTE: vapor pressure equation [Pa]
            VaPr_comp = params['vapor_pressure']

            # SECTION: kwargs
            # temperature guess [K]
            T_g0 = kwargs.get('guess_temperature', 295)

            # SECTION: activity coefficient
            # NOTE: mole fraction dict
            # convert to dict
            z_i_comp = self.__mole_fraction_comp(z_i)

            # NOTE: init model
            # init NRTL model
            activity = Activity()

            # activity inputs
            activity_inputs = kwargs.get('activity_inputs', {})

            # SECTION: optimization
            # pressure [Pa]
            P_value = P['value']

            # params
            _params = {
                'mole_fraction': z_i,
                'mole_fraction_comp': z_i_comp,
                'pressure': P_value,
                'vapor_pressure': VaPr_comp,
                'activity_model': activity_model,
                'activity': activity,
                'activity_inputs': activity_inputs
            }

            # NOTE: bubble temperature [K]
            T = None

            # check solver method
            if solver_method == 'fsolve':
                # ! fsolve
                _res0 = optimize.fsolve(
                    self.fBT, T_g0, args=(_params,), full_output=True)
                # extract
                T, infodict, ier, msg = _res0

                # check if root found
                if ier != 1:
                    raise Exception(f'root not found!, {msg}')

                # check
                if len(T) == 1:
                    T = float(T[0])
            elif solver_method == 'root':
                # ! root
                _res0 = optimize.root(
                    self.fBT, T_g0, args=(_params,))

                # check
                if not _res0.success:
                    raise Exception(f'root not found!, {_res0.message}')

                # extract
                T = _res0.x[0]

            elif solver_method == 'least-squares':
                # ! least-squares
                _res0 = optimize.least_squares(
                    self.fBT, T_g0, args=(_params,))

                # check
                if not _res0.success:
                    raise Exception(f'root not found!, {_res0.message}')

                # extract
                T = _res0.x[0]

            else:
                raise Exception('solver method not found!')

            # SECTION: vapor pressure [Pa]
            # at T (Tg)
            VaPr = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [Pa]
                eq_ = VaPr_comp[component]['value']
                args_ = VaPr_comp[component]['args']
                # update args
                args_['T'] = T
                # cal
                res_ = eq_.cal(**args_)
                # extract
                res_value_ = res_['value']
                res_unit_ = res_['unit']
                # convert to Pa
                unit_block_ = f"{res_unit_} => Pa"
                VaPr_ = pycuc.to(res_value_, unit_block_)
                # set
                VaPr[i] = VaPr_

            # SECTION: calculate activity coefficient
            # NOTE: check model
            if activity_model == 'NRTL':
                # calculate activity
                res_ = activity.NRTL(
                    self.components, z_i_comp, T, **kwargs)
                # extract
                AcCo_i = res_['value']
            elif activity_model == 'UNIQUAC':
                # calculate activity
                res_ = activity.UNIQUAC(
                    self.components, z_i_comp, T, **kwargs)
                # extract
                AcCo_i = res_['value']
            else:
                # equals unity for ideal solution
                AcCo_i = np.ones(self.component_num)

            # vapor mole fraction
            yi = np.zeros(self.component_num)
            for i in range(self.component_num):
                yi[i] = z_i[i]*AcCo_i[i]*VaPr[i]/P_value

            # NOTE: k-ratio
            K_i = np.multiply(yi, 1/z_i)

            # NOTE: results
            res = {
                "bubble_temperature": {
                    "value": float(T),
                    "unit": "K"
                },
                "pressure": {
                    "value": P_value,
                    "unit": "Pa"
                },
                "feed_mole_fraction": z_i,
                "liquid_mole_fraction": z_i,
                "vapor_mole_fraction": yi,
                "vapor_pressure": {
                    "value": VaPr,
                    "unit": "Pa"
                },
                "K_ratio": {
                    "value": K_i,
                    "unit": "dimensionless"
                },
                "activity_coefficient": {
                    "value": AcCo_i,
                    "unit": "dimensionless"
                }
            }

            # res
            return res
        except Exception as e:
            raise Exception(f"bubble temperature calculation failed! {e}")

    def fBT(self, x, params) -> float:
        '''
        bubble temperature function

        Parameters
        ----------
        x : array-like
            Guess temperature [K]
        params : dict
            Dictionary containing the following:
            - mole_fraction : liquid mole fraction (xi)
            - pressure : pressure [Pa]
            - vapor_pressure : vapor pressure equation [Pa]
        '''
        # NOTE: temperature (loop)
        T = x[0]

        # NOTE: params
        # mole fraction [array]
        z_i = params['mole_fraction']
        z_i_comp = params['mole_fraction_comp']
        # pressure [Pa]
        P = params['pressure']
        # vapor pressure calculation
        VaPr_comp = params['vapor_pressure']
        # activity model
        activity_model = params['activity_model']
        # activity
        activity: Activity = params['activity']
        # activity inputs
        activity_inputs = params['activity_inputs']

        # NOTE: check model
        if activity_model == 'NRTL':
            # calculate activity
            res_ = activity.NRTL(
                self.components, z_i_comp, T, activity_inputs=activity_inputs)
            # extract
            AcCo_i = res_['value']
        elif activity_model == 'UNIQUAC':
            # calculate activity
            res_ = activity.UNIQUAC(
                self.components, z_i_comp, T, activity_inputs=activity_inputs)
            # extract
            AcCo_i = res_['value']
        else:
            # equals unity for ideal solution
            AcCo_i = np.ones(self.component_num)

        # NOTE: calculate vapor-pressure
        # vapor pressure [Pa]
        VaPr_i = np.zeros(self.component_num)

        # looping over components
        for i, component in enumerate(self.components):
            # vapor pressure [?]
            eq_ = VaPr_comp[component]['value']
            args_ = VaPr_comp[component]['args']
            # update args
            args_['T'] = T
            # cal
            res_ = eq_.cal(**args_)
            # extract
            res_value_ = res_['value']
            res_unit_ = res_['unit']
            # convert to Pa
            unit_block_ = f"{res_unit_} => Pa"
            VaPr_ = pycuc.to(res_value_, unit_block_)

            # save
            VaPr_i[i] = VaPr_

        # NOTE: bubble pressure [Pa]
        BuPr = np.dot(z_i*AcCo_i, VaPr_i)

        # NOTE: loss
        loss = abs((P/BuPr) - 1)

        return loss

    def _DT(self, params, **kwargs):
        '''
        The `dew-point temperature` (DT) calculation determines the temperature at which the first drop of liquid condenses when a vapor mixture is cooled at a constant pressure. It identifies the temperature at which vapor will start to condense.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : vapor mole fraction (yi)
            - P: pressure [Pa]
        kwargs : dict
            additional parameters for the calculation
            - `guess_temperature`: initial guess temperature [K], default is 295 K

        Returns
        -------
        res : dict
            Dictionary containing the following:
            - dew_temperature: dew temperature [K]
            - pressure: pressure [Pa]
            - feed_mole_fraction: vapor mole fraction (yi)
            - vapor_mole_fraction: vapor mole fraction (yi)
            - liquid_mole_fraction: liquid mole fraction (xi)
            - vapor_pressure: vapor pressure [Pa]
            - activity_coefficient: activity coefficient [dimensionless]
            - K_ratio: K-ratio [dimensionless]

        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Pressure (P) and mole fraction of the components in the vapor phase (yi).
        then zi = yi.
        - Computed Information: Temperature (T) and mole fraction in the liquid phase (xi).

        The solution is obtained by the following steps:
            1. Choose an initial guess for the dew-point temperature (Tg).
            2. Calculate the vapor pressure of each component at the guessed temperature (Tg).
            3. Calculate the dew pressure (Pb) using `Raoult's law` or other methods.
            4. Check if the calculated dew pressure (Pb) equals the given pressure (P).
            5. Calculate the liquid mole fraction (xi) using `Raoult's law` or other methods.
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            z_i = params['mole_fraction']
            # pressure [Pa]
            P = params['pressure']
            # equilibrium model
            eq_model = params['equilibrium_model']
            # fugacity model
            fugacity_model = params['fugacity_model']
            # activity model
            activity_model = params['activity_model']
            # solver method
            solver_method = params['solver_method']

            # SECTION: vapor pressure calculation
            # NOTE: vapor pressure equation [Pa]
            VaPr_comp = params['vapor_pressure']

            # NOTE: kwargs
            # temperature guess [K]
            T_g0 = kwargs.get('guess_temperature', 295)

            # SECTION: activity coefficient
            # NOTE: mole fraction dict
            # convert to dict
            z_i_comp = self.__mole_fraction_comp(z_i)

            # NOTE: init model
            # init NRTL model
            activity = Activity()

            # activity inputs
            activity_inputs = kwargs.get('activity_inputs', {})

            # SECTION: optimization
            # pressure [Pa]
            P_value = P['value']

            # params
            _params = {
                'mole_fraction': z_i,
                'mole_fraction_comp': z_i_comp,
                'pressure': P_value,
                'vapor_pressure': VaPr_comp,
                'activity_model': activity_model,
                'activity': activity,
                'activity_inputs': activity_inputs
            }

            # NOTE: dew temperature [K]
            T = None

            # check solver method
            if solver_method == 'fsolve':
                # ! fsolve
                # bubble pressure [Pa]
                _res0 = optimize.fsolve(
                    self.fDT,
                    T_g0,
                    args=(_params),
                    full_output=True)

                # extract
                T, infodict, ier, msg = _res0

                # check if root found
                if ier != 1:
                    raise Exception(f'root not found!, {msg}')

                # check
                if len(T) == 1:
                    T = float(T[0])

            elif solver_method == 'root':
                # ! root
                _res0 = optimize.root(
                    self.fDT, T_g0, args=(_params,))

                # check
                if not _res0.success:
                    raise Exception(f'root not found!, {_res0.message}')

                # extract
                T = _res0.x[0]

            elif solver_method == 'least-squares':
                # ! least-squares
                _res0 = optimize.least_squares(
                    self.fDT, T_g0, args=(_params,))

                # check
                if not _res0.success:
                    raise Exception(f'root not found!, {_res0.message}')

                # extract
                T = _res0.x[0]

            else:
                raise Exception('solver method not found!')

            # NOTE: vapor pressure [Pa]
            # at T (Tg)
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [?]
                eq_ = VaPr_comp[component]['value']
                args_ = VaPr_comp[component]['args']
                # update args
                args_['T'] = T
                # cal
                res_ = eq_.cal(**args_)
                # extract
                res_value_ = res_['value']
                res_unit_ = res_['unit']
                # convert to Pa
                unit_block_ = f"{res_unit_} => Pa"
                VaPr_ = pycuc.to(res_value_, unit_block_)

                # save
                VaPr_i[i] = VaPr_

            # NOTE: calculate activity coefficient
            # NOTE: check model
            if activity_model == 'NRTL':
                # calculate activity
                res_ = activity.NRTL(
                    self.components, z_i_comp, T, **kwargs)
                # extract
                AcCo_i = res_['value']
            elif activity_model == 'UNIQUAC':
                # calculate activity
                res_ = activity.UNIQUAC(
                    self.components, z_i_comp, T, **kwargs)
                # extract
                AcCo_i = res_['value']
            else:
                # equals unity for ideal solution
                AcCo_i = np.ones(self.component_num)

            # NOTE: vapor mole fraction
            xi = np.zeros(self.component_num)
            for i in range(self.component_num):
                xi[i] = (z_i[i]*P_value)/(VaPr_i[i]*AcCo_i[i])

            # NOTE: k-ratio
            K_i = np.multiply(xi, 1/z_i)

            # NOTE: results
            res = {
                "dew_temperature": {
                    "value": float(T),
                    "unit": "K"
                },
                "pressure": {
                    "value": P_value,
                    "unit": "Pa"
                },
                "feed_mole_fraction": z_i,
                "liquid_mole_fraction": xi,
                "vapor_mole_fraction": z_i,
                "vapor_pressure": {
                    "value": VaPr_i,
                    "unit": "Pa"
                },
                "K_ratio": {
                    "value": K_i,
                    "unit": "dimensionless"
                },
                "activity_coefficient": {
                    "value": AcCo_i,
                    "unit": "dimensionless"
                }
            }

            # res
            return res
        except Exception as e:
            raise Exception(f'dew temperature calculation failed! {e}')

    def fDT(self, x, params) -> float:
        '''
        dew temperature function

        Parameters
        ----------
        x : array-like
            Guess temperature [K]
        params : dict
            Dictionary containing the following:
            - mole_fraction : liquid mole fraction (xi)
            - pressure : pressure [Pa]
            - vapor_pressure : vapor pressure equation [Pa]
        '''
        # NOTE: temperature (loop)
        T = x[0]

        # NOTE: params
        # mole fraction [array]
        z_i = params['mole_fraction']
        z_i_comp = params['mole_fraction_comp']
        # pressure [Pa]
        P = params['pressure']
        # vapor pressure calculation
        VaPr_comp = params['vapor_pressure']
        # activity model
        activity_model = params['activity_model']
        # activity
        activity: Activity = params['activity']
        # activity inputs
        activity_inputs = params['activity_inputs']

        # NOTE: calculate vapor-pressure
        # vapor pressure [Pa]
        VaPr_i = np.zeros(self.component_num)

        # looping over components
        for i, component in enumerate(self.components):
            # vapor pressure [?]
            eq_ = VaPr_comp[component]['value']
            args_ = VaPr_comp[component]['args']
            # update args
            args_['T'] = T
            # cal
            res_ = eq_.cal(**args_)
            # extract
            res_value_ = res_['value']
            res_unit_ = res_['unit']
            # convert to Pa
            unit_block_ = f"{res_unit_} => Pa"
            VaPr_ = pycuc.to(res_value_, unit_block_)

            # save
            VaPr_i[i] = VaPr_

        # NOTE: calculate activity coefficient
        # NOTE: check model
        if activity_model == 'NRTL':
            # calculate activity
            res_ = activity.NRTL(
                self.components, z_i_comp, T, activity_inputs=activity_inputs)
            # extract
            AcCo_i = res_['value']
        elif activity_model == 'UNIQUAC':
            # calculate activity
            res_ = activity.UNIQUAC(
                self.components, z_i_comp, T, activity_inputs=activity_inputs)
            # extract
            AcCo_i = res_['value']
        else:
            # equals unity for ideal solution
            AcCo_i = np.ones(self.component_num)

        # NOTE: dew pressure [Pa]
        DePr = 1/np.dot(z_i, 1/(VaPr_i*AcCo_i))

        # NOTE: loss function
        loss = abs((P/DePr) - 1)

        return loss

    def __cal_bubble_pressure(self, xi: np.ndarray, VaPr: np.ndarray) -> float:
        '''
        Calculate bubble pressure according to Raoult's law.

        Parameters
        ----------
        xi : array-like
            Liquid mole fraction of each component in the mixture.
        VaPr : array-like
            Vapor pressure of each component in the mixture.

        Returns
        -------
        BuPr : float
            Bubble pressure of the mixture.
        '''
        # bubble pressure
        BuPr = np.dot(xi, VaPr)

        # res
        return BuPr

    def __cal_dew_pressure(self, yi: np.ndarray, VaPr: np.ndarray) -> float:
        '''
        Calculate dew pressure according to Raoult's law.

        Parameters
        ----------
        yi : array-like
            Vapor mole fraction of each component in the mixture.
        VaPr : array-like
            Vapor pressure of each component in the mixture.

        Returns
        -------
        DePr : float
            Dew pressure of the mixture.
        '''
        # dew-point pressure
        DePr = 1/np.dot(yi, 1/VaPr)

        # res
        return DePr

    def _flash_checker(self,
                       z_i: List[float],
                       Pf: float,
                       Tf: float,
                       VaPr_comp: Dict) -> bool:
        '''
        Check if the flash occurs at the given pressure and temperature according to the bubble and dew pressures of the mixture; according to Raoult's law.

        Parameters
        ----------
        z_i : list
            Feed mole fraction of each component in the mixture.
        Pf : float
            Flash pressure [Pa]
        Tf : float
            Flash temperature [K]
        VaPr_comp : dict
            Dictionary containing the vapor pressure equations for each component.

        Returns
        -------
        bool
            True if the pressure and temperature are valid, False otherwise.

        Notes
        -----
        Assumption are as:

        - The system is in equilibrium.
        - The vapor pressure of each component is calculated using the provided equations.
        - The bubble pressure is calculated using Raoult's law.
        - The dew pressure is calculated using Raoult's law.
        '''
        try:
            # NOTE: vapor pressure [Pa]
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [?]
                eq_ = VaPr_comp[component]['value']
                args_ = VaPr_comp[component]['args']
                # update args
                args_['T'] = Tf
                # cal
                res_ = eq_.cal(**args_)
                # extract
                res_value_ = res_['value']
                res_unit_ = res_['unit']
                # convert to Pa
                unit_block_ = f"{res_unit_} => Pa"
                VaPr_ = pycuc.to(res_value_, unit_block_)
                # save
                VaPr_i[i] = VaPr_

            # NOTE: calculate bubble pressure
            BuPr = self.__cal_bubble_pressure(z_i, VaPr_i)
            # NOTE: calculate dew pressure
            DePr = self.__cal_dew_pressure(z_i, VaPr_i)

            # NOTE: check if the given pressure and temperature are within the valid range
            if Pf < BuPr and Pf > DePr:
                # if the pressure is between the bubble and dew pressures, return False
                return True
            else:
                # if the pressure is outside the valid range, return False
                return False
        except Exception as e:
            # if any error occurs, return False
            raise Exception(f"Error in flash checker: {e}")

    def _IFL(self, params, **kwargs):
        '''
        The `isothermal-flash` (IFL) calculation This calculation determines the vapor and liquid phase compositions and amounts at a specified temperature and pressure. The system is "flashed" isothermally, meaning the temperature is kept constant while the phase behavior is calculated for a mixture.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : feed mole fraction (zi)
            - P: pressure [Pa]
            - T: temperature [K]
        kwargs : dict
            additional parameters for the calculation

        Returns
        -------
        res : dict
            Dictionary containing the following:
            - vapor_to_liquid_ratio: V/F ratio [dimensionless]
            - liquid_to_vapor_ratio: L/F ratio [dimensionless]
            - feed_mole_fraction: liquid mole fraction (zi)
            - liquid_mole_fraction: liquid mole fraction (xi)
            - vapor_mole_fraction: vapor mole fraction (yi)
            - vapor_pressure: vapor pressure [Pa]
            - K_ratio: K-ratio [dimensionless]
            - temperature: temperature [K]
            - pressure: pressure [Pa]

        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Temperature (T), pressure (P), and mole fraction of the components in the liquid phase (zi).
        - Computed Information: Vapor-to-liquid ratio (V/F), liquid mole fraction (xi), vapor mole fraction (yi), and liquid-to-vapor ratio (L/F).

        The solution is obtained by the following steps:
            1. Calculate the K ratio (Ki) using Raoult's law.
            2. Choose an initial guess for the vapor-to-liquid ratio (V/F).
            3. Solve a system of non-linear equations to find the V/F ratio that satisfies the mass balance equations.
            4. Calculate the liquid and vapor mole fractions (xi and yi) using the V/F ratio and K ratio.
            5. Calculate the liquid-to-vapor ratio (L/F) from the V/F ratio.
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            z_i = params['mole_fraction']
            # pressure [Pa]
            P = params['pressure']
            # temperature [K]
            T = params['temperature']
            # equilibrium model
            eq_model = params['equilibrium_model']
            # fugacity model
            fugacity_model = params['fugacity_model']
            # activity model
            activity_model = params['activity_model']
            # solver method
            solver_method = params['solver_method']

            # SECTION: vapor pressure calculation
            # NOTE: vapor pressure equation [Pa]
            VaPr_comp = params['vapor_pressure']

            # NOTE: kwarg

            # NOTE: set values
            # pressure [Pa]
            P_value = P['value']
            # temperature [K]
            T_value = T['value']

            # NOTE: ki ratio (Raoult's law)
            # pressure and temperature are constant (Raoutl's law)
            K_i = np.zeros(self.component_num)

            # vapor pressure [Pa]
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [?]
                eq_ = VaPr_comp[component]['value']
                args_ = VaPr_comp[component]['args']
                # update args
                args_['T'] = T_value
                # cal
                res_ = eq_.cal(**args_)
                # extract
                res_value_ = res_['value']
                res_unit_ = res_['unit']
                # convert to Pa
                unit_block_ = f"{res_unit_} => Pa"
                VaPr_ = pycuc.to(res_value_, unit_block_)
                # set
                VaPr_i[i] = VaPr_
                K_i[i] = VaPr_/P_value

            # SECTION: optimization
            # NOTE: set params
            _params = {
                'mole_fraction': z_i,
                'K_ratio': K_i
            }

            # NOTE: check solver method
            if solver_method == 'least_squares':
                # ! least-squares
                # NOTE: initial guess
                N = self.component_num
                # Initial guess: beta, x1, x2, ..., x_{N-1}
                x0 = [0.5] + [1.0 / N] * (N - 1)

                # NOTE: Bounds
                # Small value to avoid numerical issues with log(0), etc.
                eps = kwargs.get('eps', 1e-8)

                lower_bounds = [eps] + [eps] * \
                    (N - 1)         # β and each x_i ≥ eps
                upper_bounds = [1.0 - eps] + [1.0] * \
                    (N - 1)   # β ≤ 1−eps, x_i ≤ 1

                bounds = (lower_bounds, upper_bounds)

                # V/F
                _res = optimize.least_squares(
                    self.fIFL, x0, args=(_params,), bounds=bounds)

                # check if root found
                if _res.success is False:
                    raise Exception(f'root not found!, {_res.message}')

                # -> result analysis
                V_F_ratio = _res.x[0]
                x_liq = np.zeros_like(z_i)
                x_liq[:-1] = _res.x[1:]
                x_liq[-1] = 1 - np.sum(x_liq[:-1])

            elif solver_method == 'minimize':
                # ! minimize
                # NOTE: initial guess
                N = self.component_num
                # Initial guess: beta, x1, x2, ..., x_{N-1}
                x0 = [0.5] + [1.0 / N] * (N - 1)

                # NOTE: Bounds
                # Small value to avoid numerical issues with log(0), etc.
                eps = kwargs.get('eps', 1e-8)

                # individual bounds (min, max)
                bounds = [(eps, 1.0 - eps)] + [(eps, 1.0)
                                               for _ in range(N - 1)]

                # V/F
                _res = optimize.minimize(
                    self.fIFL2,
                    x0=x0,
                    args=(_params,),
                    constraints=self.flash_constraints(z_i, K_i),
                    bounds=bounds,
                )

                # check if root found
                if _res.success is False:
                    raise Exception(f'root not found! {_res.message}')

                # -> result analysis
                V_F_ratio = _res.x[0]
                x_liq = np.zeros_like(z_i)
                x_liq[:-1] = _res.x[1:]
                x_liq[-1] = 1 - np.sum(x_liq[:-1])

            else:
                raise Exception('solver method not found!')

            # SECTION: liquid/vapor mole fraction
            xy_ = self.xy_flash(V_F_ratio, z_i, K_i)
            # set
            xi = xy_['liquid']
            yi = xy_['vapor']

            # NOTE: calculate L/F
            L_F_ratio = 1 - V_F_ratio

            # NOTE: results
            res = {
                "V_F_ratio": {
                    "value": float(V_F_ratio),
                    "unit": "dimensionless"
                },
                "L_F_ratio": {
                    "value": float(L_F_ratio),
                    "unit": "dimensionless"
                },
                "feed_mole_fraction": z_i,
                "liquid_mole_fraction": xi,
                "vapor_mole_fraction": yi,
                "vapor_pressure": {
                    "value": VaPr_i,
                    "unit": "Pa"
                },
                "K_ratio": {
                    "value": K_i,
                    "unit": "dimensionless"
                },
                "temperature": {
                    "value": T_value,
                    "unit": "K"
                },
                "pressure": {
                    "value": P_value,
                    "unit": "Pa"
                },
            }

            # res
            return res
        except Exception as e:
            raise Exception(f"flash isothermal failed! {e}")

    def fIFL(self, x, params):
        '''
        Flash isothermal function, according to Raoult's law.

        Parameters
        ----------
        x : array-like
            VF : vapor-feed ratio [dimensionless]
            x_i : liquid mole fraction [dimensionless]
        params : tuple
            Tuple containing the following:
            - zi : feed mole fraction [dimensionless]
            - Ki : K ratio [dimensionless]
        '''
        # NOTE: extract variables
        # vapor fraction
        VF = x[0]
        # liquid mole fraction
        x_i = np.zeros(self.component_num)
        x_i[:-1] = x[1:]
        # last x_i
        x_i[-1] = 1 - np.sum(x_i[:-1])

        # NOTE: params
        # feed mole fraction
        z_i = params['mole_fraction']
        # K ratio
        K_i = params['K_ratio']

        # NOTE: vapor mole fraction
        y_i = K_i * x_i

        # SECTION: check optimization region
        # Avoid division by zero or invalid regions
        if VF <= 0 or VF >= 1:
            return 1e6  # Return a large value to indicate invalid region

        # NOTE: function
        eqs = []

        # looping over components
        for i in range(self.component_num - 1):
            eq = z_i[i] - ((1 - VF) * x_i[i] + VF * y_i[i])
            eqs.append(eq)

        # Closure relations
        eqs.append(np.sum(x_i) - 1)   # sum x_i = 1
        eqs.append(np.sum(y_i) - 1)   # sum y_i = 1

        return eqs

    def xy_flash(self, V_F_ratio: float, z_i: np.ndarray, K_i: np.ndarray) -> Dict[str, np.ndarray]:
        '''
        Calculate liquid/vapor mole fraction (xi, yi) using V/F ratio and K ratio.

        Parameters
        ----------
        V_F_ratio : float
            Vapor-to-liquid ratio (V/F).
        zi : array-like
            Feed mole fraction of each component in the mixture.
        Ki : array-like
            K ratio of each component in the mixture.

        Returns
        -------
        dict
            Dictionary containing the following:
            - liquid: liquid mole fraction (xi) [dimensionless]
            - vapor: vapor mole fraction (yi) [dimensionless]
        '''
        try:
            # init
            x_i = np.zeros(self.component_num)
            y_i = np.zeros(self.component_num)

            # liquid/vapor mole fraction
            for i in range(self.component_num):
                x_i[i] = (z_i[i])/(1+(V_F_ratio)*(K_i[i]-1))
                y_i[i] = K_i[i]*x_i[i]

            # res
            return {
                'liquid': x_i,
                'vapor': y_i
            }
        except Exception as e:
            raise Exception(f"Error in xy_flash calculation: {e}")

    def fIFL2(self, x, params):
        '''
        Flash isothermal function (nonlinear equations), according to Raoult's law.

        Parameters
        ----------
        x : array-like
            VF : vapor-feed ratio [dimensionless]
            x_i : liquid mole fraction [dimensionless]
        params : tuple
            Tuple containing the following:
            - zi: feed mole fraction (zi)
            - Ki: K ratio of each component in the mixture
        '''
        # NOTE: extract variables
        # vapor fraction
        VF = x[0]
        # liquid mole fraction
        x_i = np.zeros(self.component_num)
        x_i[:-1] = x[1:]
        # last x_i
        x_i[-1] = 1 - np.sum(x_i[:-1])

        # NOTE: params
        # feed mole fraction
        z_i = params['mole_fraction']
        # K ratio
        K_i = params['K_ratio']

        # NOTE: vapor mole fraction
        y_i = K_i * x_i

        # SECTION: system of nonlinear equations (NLE)
        # Residuals of component balances
        residuals = z_i - ((1 - VF) * x_i + VF * y_i)

        # Objective: minimize sum of squared residuals
        return np.sum(residuals**2)

    def flash_constraints(self, z, K):
        N = len(z)

        def liquid_sum(vars):
            x = np.zeros(N)
            x[:-1] = vars[1:]
            x[-1] = 1.0 - np.sum(x[:-1])
            return np.sum(x) - 1

        def vapor_sum(vars):
            x = np.zeros(N)
            x[:-1] = vars[1:]
            x[-1] = 1.0 - np.sum(x[:-1])
            y = K * x
            return np.sum(y) - 1

        return [
            {'type': 'eq', 'fun': liquid_sum},
            {'type': 'eq', 'fun': vapor_sum},
        ]

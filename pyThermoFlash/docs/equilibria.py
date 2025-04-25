# VAPOR-LIQUID EQUILIBRIUM (VLE) CALCULATIONS
# --------------------------------------------
# import libs
from typing import List, Dict, Optional
import numpy as np
from scipy import optimize
# local


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

    def BP(self, params, **kwargs):
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

        Returns
        -------
        dict
            Dictionary containing the following:
            - bubble_pressure: bubble pressure [Pa]
            - temperature: temperature [K]
            - feed_mole_fraction: liquid mole fraction [dimensionless]
            - vapor_mole_fraction: vapor mole fraction [dimensionless]
            - vapor_pressure: vapor pressure [Pa]
            - activity_coefficient: activity coefficient [dimensionless]
            - K_ratio: K ratio [dimensionless]

        Notes
        -----
        The summary of the calculation is as follows:

        - Known Information: Temperature (T) and mole fraction of the components in the liquid phase (xi).
        then zi = xi.
        - Computed Information: Pressure (P) and mole fraction in the vapor phase (yi).

        - The solution is obtained by the following steps:
            1. Calculate the vapor pressure of each component at the given temperature (T).
            2. Calculate the bubble pressure (P) using Raoult's law.
            3. Calculate the mole fraction in the vapor phase (yi) using `Raoult's law`.
        '''
        try:
            # SECTION: data
            # NOTE: params
            # mole fraction [array]
            z_i = params['mole_fraction']
            # temperature [K]
            T = params['temperature']
            # activity model
            activity_model = params.get('activity_model')

            # NOTE: vapor pressure [Pa]
            VaPr_comp = params['vapor_pressure']

            # SECTION: activity coefficient
            if activity_model == 'NRTL':
                # init NRTL model
                AcCo_i = np.zeros(self.comp_num)
            elif activity_model == 'UNIQUAC':
                # init UNIQUAC model
                AcCo_i = np.zeros(self.comp_num)
            else:
                # equals unity for ideal solution
                AcCo_i = np.ones(self.comp_num)

            # SECTION: vapor pressure calculation
            # vapor pressure [Pa]
            VaPr_i = np.zeros(self.comp_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [Pa]
                VaPr_i[i] = VaPr_comp[component]

            # bubble pressure [Pa]
            BuPr = np.sum(AcCo_i*z_i*VaPr_i)

            # vapor mole fraction
            y_i = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                y_i[i] = z_i[i]*VaPr_i[i]*AcCo_i[i]/BuPr

            # NOTE: Ki ratio
            K_i = np.multiply(y_i, 1/z_i)

            # NOTE: results
            res = {
                "bubble_pressure": {
                    "value": BuPr,
                    "unit": "Pa"
                },
                "temperature": {
                    "value": T,
                    "unit": "K"
                },
                "feed_mole_fraction": z_i,
                "vapor_mole_fraction": y_i,
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

    def DP(self, params, **kwargs):
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
            # activity model
            activity_model = params.get('activity_model')

            # NOTE: vapor pressure [Pa]
            VaPr_comp = params['vapor_pressure']

            # vapor pressure [Pa]
            VaPr_i = np.zeros(self.component_num)

            # looping over components
            for i, component in enumerate(self.components):
                # vapor pressure [Pa]
                VaPr_i[i] = VaPr_comp[component]

            # NOTE: dew pressure [Pa]
            DePr = 1/np.dot(y_i, 1/VaPr_i)

            # NOTE: liquid mole fraction
            x_i = np.zeros(self.component_num)

            # looping over components
            for i in range(self.component_num):
                # mole fraction
                x_i[i] = y_i[i]*DePr/VaPr_i[i]

            # NOTE: k-ratio
            K_i = np.multiply(y_i, 1/x_i)

            # res
            res = {
                "dew_pressure": {
                    "value": DePr,
                    "unit": "Pa"
                },
                "temperature": {
                    "value": T,
                    "unit": "K"
                },
                "vapor_mole_fraction": y_i,
                "liquid_mole_fraction": x_i,
                "vapor_pressure": {
                    "value": VaPr_i,
                    "unit": "Pa"
                },
                "activity_coefficient": {
                    "value": np.ones(self.comp_num),
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

    def BT(self, params, **kwargs):
        '''
        The bubble-point temperature (BT) calculation determines the temperature at which the first bubble of vapor forms when a liquid mixture is heated at a constant pressure. It helps identify the temperature at which the liquid will start vaporizing.

        Parameters
        ----------
        params : dict
            Dictionary containing the following:
            - zi : liquid mole fraction (xi)
            - P: pressure [Pa]
        kwargs : dict
            additional parameters for the calculation

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
            zi = params['mole_fraction']
            # pressure [Pa]
            P = params['pressure']
            # activity model
            activity_model = params['activity_model']

            # NOTE: config
            Tg0 = kwargs.get('Tg0', 295)

            # params
            _params = (self.comp_num, zi, P, VaPeCal)

            # NOTE: bubble temperature [K]
            _res0 = optimize.fsolve(self.btFunction, Tg0, args=(_params,))
            # ->
            # REVIEW
            T = _res0[0]

            # NOTE: vapor pressure [Pa]
            # at T (Tg)
            VaPe = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                # REVIEW
                VaPe[i] = self.pool[i].vapor_pressure(T, VaPeCal)

            # vapor mole fraction
            yi = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                yi[i] = zi[i]*VaPe[i]/P

            _res = {
                "T": T,
                "yi": yi,
                'xi': zi,
                "VaPe": VaPe
            }

            # res
            return _res
        except Exception as e:
            raise Exception(e)

    def btFunction(self, x, params):
        '''
        args:
            x: guess temperature *** array *** [K]
        '''
        # Tg
        Tg = x[0]
        # params
        compNo, zi, P, VaPeCal = params
        # calculate vapor-pressure

        # vapor pressure [Pa]
        VaPe = np.zeros(compNo)
        for i in range(compNo):
            # REVIEW
            VaPe[i] = self.pool[i].vapor_pressure(Tg, VaPeCal)

        # bubble pressure [Pa]
        BuPr = np.dot(zi, VaPe)

        # loss
        loss = abs((P/BuPr) - 1)

        return loss

    def cal_dew_temperature(self, params, config):
        '''
        dew temperature calculation

        args:
            params:
                1. feed mole fraction *** array *** (zi=xi)
                2. system pressure [Pa]
            config:
                1. Tg0: initial guess temperature
                2. VaPeCal: vapor-pressure calculation method (default: polynomial)

        knowns:
            1. P
            2. z[i] = y[i]

        unknowns:
            1. T (dew temperature)
            2. y[i]

        solutions:
            1. guess Tg
            2. cal P[i,sat] at Tg
            3. cal bubble pressure Pb at Tg
            4. check Pb equals P
            4. cal x[i]
        '''
        try:
            # params
            zi = params.get('zi', [])
            P = params.get('P', 0)

            # config
            VaPeCal = config.get('VaPeCal', 'polynomial')
            Tg0 = config.get('Tg0', 295)

            # params
            _params = (self.comp_num, zi, P, VaPeCal)
            # bubble pressure [Pa]
            _res0 = optimize.fsolve(self.dtFunction, Tg0, args=(_params,))
            # ->
            # REVIEW
            Tg = _res0[0]

            # vapor pressure [Pa]
            # at T (Tg)
            VaPe = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                # REVIEW
                VaPe[i] = self.pool[i].vapor_pressure(Tg, VaPeCal)

            # vapor mole fraction
            xi = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                xi[i] = zi[i]*P/VaPe[i]

            # res
            return Tg, xi, VaPe
        except Exception as e:
            raise Exception(e)

    def dtFunction(self, Tg, params):
        '''
        args:
            Tg: guess temperature [K]
        '''
        # params
        compNo, zi, P, VaPeCal = params
        # calculate vapor-pressure

        # vapor pressure [Pa]
        VaPe = np.zeros(compNo)
        for i in range(compNo):
            # REVIEW
            VaPe[i] = self.pool[i].vapor_pressure(Tg, VaPeCal)

        # dew pressure [Pa]
        DePr = np.dot(zi, VaPe)

        # loss
        loss = abs((P/DePr) - 1)

        return loss

    def vapor_pressure_mixture(self, T, mode):
        '''
        calculate mixture vapor-pressure
        '''
        # vapor pressure [Pa]
        VaPe = np.zeros(self.comp_num)
        for i in range(self.comp_num):
            # REVIEW
            VaPe[i] = self.pool[i].vapor_pressure(T, mode)

        # res
        return VaPe

    def bubble_pressure(self, xi, VaPe):
        '''
        calculate bubble pressure

        args:
            xi: liquid mole fraction
            VaPe: vapor-pressure [Pa]
        '''
        # bubble pressure
        BuPr = np.dot(xi, VaPe)

        # res
        return BuPr

    def dew_pressure(self, yi, VaPe):
        '''
        calculate dew pressure
        '''
        # dew-point pressure
        DePr = 1/np.dot(yi, 1/VaPe)

        # res
        return DePr

    def cal_flash_isothermal(self, params, config):
        '''
        isothermal flash calculation

        knowns:
            1. zi
            2. P
            3. T

        unknowns:
            1. xi
            2. yi
            3. V
            4. L
        '''
        try:
            # params
            zi = params.get('zi', [])
            P_flash = params.get('P_flash', 0)
            T_flash = params.get('T_flash', 0)
            VaPri = params.get('VaPe', [])

            # config
            VaPeCal = config.get('VaPeCal', 'polynomial')
            V_F_ratio_g0 = config.get('guess_V_F_ratio', 0.5)

            # ki ratio (Raoult's law)
            Ki = VaPri/P_flash

            # params
            _params = (self.comp_num, zi, Ki)
            # V/F
            _res0 = optimize.fsolve(
                self.fitFunction, V_F_ratio_g0, args=(_params,))
            # ->
            V_F_ratio = _res0[0]

            # liquid/vapor mole fraction
            xi, yi = self.xyFlash(self.comp_num, V_F_ratio, zi, Ki)

            # L/F
            L_F_ratio = 1 - V_F_ratio

            # res
            return V_F_ratio, L_F_ratio, xi, yi

        except Exception as e:
            raise Exception("flash isothermal failed!")

    def fitFunction(self, x, params):
        '''
        flash isothermal function

        args:
            x: V/F guess
            params:
                zi: feed mole fraction
                Ki: K ratio
        '''
        # V/F
        V_F_ratio = x[0]

        # params
        compNo, zi, Ki = params

        fi = np.zeros(compNo)
        for i in range(compNo):
            fi[i] = (zi[i]*(1-Ki[i]))/(1+(V_F_ratio)*(Ki[i]-1))

        f = np.sum(fi)

        return f

    def xyFlash(self, compNo, V_F_ratio, zi, Ki):
        '''
        calculate liquid/vapor mole fraction
        '''
        xi = np.zeros(compNo)
        yi = np.zeros(compNo)

        for i in range(compNo):
            xi[i] = (zi[i])/(1+(V_F_ratio)*(Ki[i]-1))
            yi[i] = Ki[i]*xi[i]

        # res
        return xi, yi

    def flashIsothermalV2(self, params, config):
        '''
        isothermal flash calculation

        knowns:
            1. zi
            2. P
            3. T

        unknowns:
            1. xi
            2. yi
            3. V
            4. L
        '''
        try:
            # params
            F = params.get('F', 0)
            zi = params.get('zi', [])
            P_flash = params.get('P_flash', 0)
            T_flash = params.get('T_flash', 0)
            VaPri = params.get('VaPe', [])

            # config
            VaPeCal = config.get('VaPeCal', 'polynomial')

            # ki ratio (Raoult's law)
            Ki = VaPri/P_flash

            # unknown no
            unknownNo = self.comp_num + 2

            # initial guess
            L0 = 0.001
            V0 = 0.001
            xi0 = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                xi0[i] = 0.01
            # set
            _var0 = [L0, V0, *xi0]

            # bounds
            bU = []
            bL = []
            bounds = []
            # lower
            for i in range(unknownNo):
                _bl = 0
                bL.append(_bl)

            # upper
            for i in range(unknownNo):
                _bu = 0.99
                bU.append(_bu)

            bounds = [bL, bU]

            # params
            _params = (self.comp_num, zi, F, Ki)
            # a system of non-linear equation
            # _res0 = optimize.fsolve(
            #     self.fitSystemFunction, _var0, args=(_params,))

            _res0 = optimize.least_squares(
                self.fitSystemFunction, _var0, args=(_params,), bounds=bounds)

            # ! check
            if _res0.success is False:
                raise Exception('root not found!')

            # sol
            x = _res0.x
            # ->
            L_sol = x[0]
            V_sol = x[1]
            V_F_ratio = V_sol/F

            # liquid/vapor mole fraction
            xi, yi = self.xyFlash(self.comp_num, V_F_ratio, zi, Ki)

            # L/F
            L_F_ratio = 1 - V_F_ratio

            # res
            return V_F_ratio, L_F_ratio, xi, yi

        except Exception as e:
            raise Exception("flash isothermal failed!")

    def fitSystemFunction(self, x, params):
        '''
        flash isothermal function (nonlinear equations)

        unknowns:
            1. L
            2. V
            3. x[i]

        args:
            x:
                1. L
                2. V
                3. x[i]
            params:
                compNo: component number
                zi: feed mole fraction
                F: feed flowrate

        '''
        # params
        compNo, zi, F, Ki = params

        # xi
        xi = np.zeros(compNo)

        # set x
        L = x[0]
        V = x[1]
        for i in range(compNo):
            xi[i] = x[i+2]

        # system of nonlinear equations (NLE)
        fi = np.zeros(compNo+2)
        # overall mole balance
        fi[0] = F - (L + V)
        # component mole balance (without reaction)
        for i in range(compNo):
            # solve the second NLE
            fi[i+1] = F*zi[i] - (L*xi[i] + V*xi[i]*Ki[i])
        # constraint
        fi[-1] = np.sum(xi) - 1

        return fi

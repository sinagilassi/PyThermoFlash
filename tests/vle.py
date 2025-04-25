# VAPOR-LIQUID EQUILIBRIUM (VLE) CALCULATIONS
# --------------------------------------------
# import libs
from math import exp, log
import numpy as np
from scipy import optimize
# local
from ..config import (
    MODIFIED_RAOULT_MODEL, R_CONST, VAN_LAAR_ACTIVITY_MODEL, WILSON_ACTIVITY_MODEL
)
from PyCTPM.docs.dThermo import ModifiedRackettEquation
from PyCTPM.docs.excessproperties import ExcessProperties
from PyCTPM.docs.eos import eosClass
from PyCTPM.docs.margules import Margules
from PyCTPM.docs.activity import ActivityClass


class VLEClass(ExcessProperties, Margules, ActivityClass):
    '''
    vapor-liquid equilibrium calculation
    '''

    def __init__(self, pool):
        self.pool = pool
        # set
        self.compNo = len(pool)

        # init class
        # ExcessProperties.__init__(self, pool)
        # Margules.__init__(self, pool)
        # ActivityClass.__init__(self, pool)

    def cal_bubble_pressure(self, params, config):
        '''
        Bubble pressure calculation

        Parameters
        ----------
        params : dict
            1. zi: liquid mole fraction
            2. T: temperature [K]
        config : dict
            1. VaPeCal: vapor pressure calculation method
            2. model: vle model
            3. AcCoModel: activity coefficient model
                1. name:
                    a) van laar: using eos to predict activity coefficient
                    b) wilson: needs parameters
                2. params:
                    a) wilson: aij ~ alpha_ij


        Returns
        -------



        Notes
        -----
        knowns:
            1. T
            2. z[i] = x[i]

        unknowns:
            1. P
            2. y[i]

        solutions:
            1. cal P[i,sat] at T
            2. cal bubble pressure
            3. cal y[i]
        '''
        try:
            # params
            zi = params.get('zi', [])
            T = params.get('T', 0)

            # config
            VaPeCal = config.get('VaPeCal')
            model = config.get('model')
            AcCoModel = config.get('AcCoModel')
            AcCoModelName = AcCoModel.get('name', 0)
            AcCoModelParameters = AcCoModel.get('params', 0)

            # activity coefficient
            if model == MODIFIED_RAOULT_MODEL:
                # !check
                if AcCoModelName == VAN_LAAR_ACTIVITY_MODEL:
                    # ! using eos to calculate activity-model
                    # set eos class
                    ai, bi = eosClass.abVDM(self.pool)
                    # set excess properties class
                    # ExcessPropertiesClass = ExcessProperties(self.pool)
                    # activity coefficient
                    AcCo = self.VanLaar_activity_coefficient(
                        zi, ai, bi, T)
                elif AcCoModelName == WILSON_ACTIVITY_MODEL:
                    AcCo = self.Wilson_activity_coefficient(
                        zi, T, AcCoModelParameters)
            else:
                # equals unity for ideal solution
                AcCo = np.ones(self.compNo)

            # vapor pressure [Pa]
            VaPe = np.zeros(self.compNo)
            for i in range(self.compNo):
                # REVIEW
                VaPe[i] = self.pool[i].vapor_pressure(T, mode=VaPeCal)

            # bubble pressure [Pa]
            BuPr = np.sum(AcCo*zi*VaPe)

            # vapor mole fraction
            yi = np.zeros(self.compNo)
            for i in range(self.compNo):
                yi[i] = zi[i]*VaPe[i]*AcCo[i]/BuPr

            # Ki ratio
            Ki = np.multiply(yi, 1/zi)

            res = {
                "P": BuPr,
                "T": T,
                "yi": yi,
                'xi': zi,
                "VaPe": VaPe,
                "AcCo": AcCo
            }

            # res
            return res
        except Exception as e:
            raise Exception(e)

    def dewPressure(self, params, config):
        '''
        dew pressure calculation

        knowns:
            1. T
            2. z[i] = y[i]

        unknowns:
            1. P
            2. x[i]

        solutions:
            1. cal P[i,sat] at T
            2. cal bubble pressure
            3. cal x[i]
        '''
        try:
            # params
            T, zis = params

            # config
            VaPeCal = config['VaPeCal']

            # compo no
            compNo = zis.size

            # vapor pressure [Pa]
            VaPe = np.zeros(compNo)
            for i in range(compNo):
                # REVIEW
                VaPe[i] = 1

            # dew pressure [Pa]
            DePr = 1/np.dot(zis, 1/VaPe)

            # liquid mole fraction
            xis = np.zeros(compNo)
            for i in range(compNo):
                xis[i] = zis[i]*DePr/VaPe[i]

            # res
            return xis, DePr
        except Exception as e:
            raise Exception(e)

    def bubbleTemperature(self, params, config):
        '''
        bubble temperature calculation

        args:
            params:
                1. feed mole fraction *** array *** (zi=xi)
                2. system pressure [Pa]
            config:
                1. Tg0: initial guess temperature
                2. VaPeCal: vapor-pressure calculation method (default: polynomial)

        knowns:
            1. P
            2. z[i] = x[i]

        unknowns:
            1. T
            2. y[i]

        solutions:
            1. guess Tg
            2. cal P[i,sat] at Tg
            3. cal bubble pressure Pb at Tg
            4. check Pb equals P
            4. cal y[i]
        '''
        try:
            # params
            zi = params.get('zi', [])
            P = params.get('P', 0)

            # config
            VaPeCal = config.get('VaPeCal', 'polynomial')
            Tg0 = config.get('Tg0', 295)

            # params
            _params = (self.compNo, zi, P, VaPeCal)
            # bubble temperature [K]
            _res0 = optimize.fsolve(self.btFunction, Tg0, args=(_params,))
            # ->
            # REVIEW
            T = _res0[0]

            # vapor pressure [Pa]
            # at T (Tg)
            VaPe = np.zeros(self.compNo)
            for i in range(self.compNo):
                # REVIEW
                VaPe[i] = self.pool[i].vapor_pressure(T, VaPeCal)

            # vapor mole fraction
            yi = np.zeros(self.compNo)
            for i in range(self.compNo):
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

    def dewTemperature(self, params, config):
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
            _params = (self.compNo, zi, P, VaPeCal)
            # bubble pressure [Pa]
            _res0 = optimize.fsolve(self.dtFunction, Tg0, args=(_params,))
            # ->
            # REVIEW
            Tg = _res0[0]

            # vapor pressure [Pa]
            # at T (Tg)
            VaPe = np.zeros(self.compNo)
            for i in range(self.compNo):
                # REVIEW
                VaPe[i] = self.pool[i].vapor_pressure(Tg, VaPeCal)

            # vapor mole fraction
            xi = np.zeros(self.compNo)
            for i in range(self.compNo):
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

    def vaporPressureMixture(self, T, mode):
        '''
        calculate mixture vapor-pressure
        '''
        # vapor pressure [Pa]
        VaPe = np.zeros(self.compNo)
        for i in range(self.compNo):
            # REVIEW
            VaPe[i] = self.pool[i].vapor_pressure(T, mode)

        # res
        return VaPe

    def calBubblePressure(self, xi, VaPe):
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

    def calDewPressure(self, yi, VaPe):
        '''
        calculate dew pressure
        '''
        # dew-point pressure
        DePr = 1/np.dot(yi, 1/VaPe)

        # res
        return DePr

    def flashIsothermal(self, params, config):
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
            _params = (self.compNo, zi, Ki)
            # V/F
            _res0 = optimize.fsolve(
                self.fitFunction, V_F_ratio_g0, args=(_params,))
            # ->
            V_F_ratio = _res0[0]

            # liquid/vapor mole fraction
            xi, yi = self.xyFlash(self.compNo, V_F_ratio, zi, Ki)

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
            unknownNo = self.compNo + 2

            # initial guess
            L0 = 0.001
            V0 = 0.001
            xi0 = np.zeros(self.compNo)
            for i in range(self.compNo):
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
            _params = (self.compNo, zi, F, Ki)
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
            xi, yi = self.xyFlash(self.compNo, V_F_ratio, zi, Ki)

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

    def activityCoefficientUsingModifiedRaoult(self, xi, yi, P, VaPr):
        '''
        calculate activity coefficient using modified raoult's law

        args:
            xi: liquid mole fraction
            yi: vapor mole fraction
            P: pressure [Pa]
            VaPr: vapor pressure [Pa]
        '''
        # AcCo
        AcCo = np.zeros(self.compNo)
        for i in range(self.compNo):
            AcCo[i] = (yi[i]*P)/(xi[i]*VaPr[i])

        # res
        return AcCo

    def MargulesParameterObjectiveFunction(self, x, params):
        '''
        Margules 1-parameter function

            args:
                x: Aij *** array ***
                params:
        '''
        # params
        xi_exp, ExMoGiEn_exp, parameterNo = params

        # number of experimental data
        expDataNo = xi_exp.shape[0]

        # calculate excess molar gibbs energy
        ExMoGiEn_cal = np.zeros(expDataNo)
        for i in range(expDataNo):
            # xi
            _xi = xi_exp[i, :]

            # calculate activity coefficient
            _AcCo = self.Margules_activity_coefficient(_xi, x)

            ExMoGiEn_cal[i] = self.ExcessMolarGibbsFreeEnergy(_xi, _AcCo)

        # obj function
        return ExMoGiEn_exp - ExMoGiEn_cal

    def margulesParameterEstimator(self, params):
        '''
        parameter estimation of Margules model for a binary system

        args:
            params:
                1. liquid mole fraction
                2. excess molar gibbs energy
                3. margules model (1-2 parameter)

        '''
        # ! check
        xi_exp, ExMoGiEn, parameterNo = params

        # initial guess
        if parameterNo == 1:
            A0 = 0.5
        elif parameterNo == 2:
            A0 = [0.5, 0.5]
        else:
            A0 = 0
            raise Exception("check A0")

        res = optimize.least_squares(
            self.MargulesParameterObjectiveFunction, A0, args=(params,))

        return res

    def WilsonParameterObjectiveFunction(self, x, params):
        '''
        Wilson function

            args:
                x: Aij *** array ***
                params:
        '''
        # params
        xi_exp, ExMoGiEn_exp = params

        # number of experimental data
        expDataNo = xi_exp.shape[0]

        # calculate excess molar gibbs energy
        ExMoGiEn_cal = np.zeros(expDataNo)
        for i in range(expDataNo):
            # xi
            _xi = xi_exp[i, :]

            # calculate activity coefficient
            _AcCo = self.wilson_activity_coefficient_parameter_estimation(
                _xi, x)

            ExMoGiEn_cal[i] = self.ExcessMolarGibbsFreeEnergy(_xi, _AcCo)

        # obj function
        return ExMoGiEn_exp - ExMoGiEn_cal

    def WilsonTemperatureIndependentParametersFunction(self, x, params):
        '''
        find temperature-independent parameters (alpha)

        *** alpha is constant ***

        args:
            x: alpha[i,j]
            params:
                1. MoVoi: molar-volume [m^3/mol]
                2. Aij: temperature-dependent parameters
                3. T: temperature [K]
        '''
        # params
        MoVoi, Aij, T = params

        # vars
        alpha_ij = np.zeros((self.compNo, self.compNo))
        y = []
        k = 0

        # set
        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    alpha_ij[i, j] = x[k]
                    k += 1

        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    _y = Aij[i, j] - (MoVoi[j]/MoVoi[i]) * \
                        exp(-1*alpha_ij[i, j]/(R_CONST*T))
                    y.append(_y)

        # res
        return y

    def WilsonParameterEstimator(self, params):
        '''
        parameter estimation of wilson model for a binary system

        the parameters are:
            1. Aij (All cross parameters are equal to each other)
                unknownNo: component number
            2. alpha_ij (is temperature independent)

        args:
            params:
                1. liquid mole fraction
                2. excess molar gibbs energy

        '''
        # ! check
        xi_exp, ExMoGiEn, T = params

        # vars
        MoVoi = np.zeros(self.compNo)

        # molar volume [m^3/mol]
        for i in range(self.compNo):
            _comp = self.components[i]
            _Pc = float(_comp.Pc)*1e5
            _Tc = _comp.Tc
            _w = _comp.w
            MoVoi[i] = ModifiedRackettEquation(T, _Pc, _Tc, _w)

        #! least-square function for Aij
        fun1 = self.WilsonParameterObjectiveFunction

        # initial guess (number of unknown parameters)
        A0 = 0.1*np.ones(self.compNo)

        # params
        params1 = (xi_exp, ExMoGiEn)

        # bounds
        unknownNo = self.compNo
        bU = []
        bL = []
        bounds = []
        # lower
        for i in range(unknownNo):
            _bl = 0
            bL.append(_bl)

        # upper
        for i in range(unknownNo):
            _bu = 3
            bU.append(_bu)

        bounds = [bL, bU]

        res0 = optimize.least_squares(
            fun1, A0, args=(params1,), bounds=bounds)

        # temperature-dependent parameters
        Aij = np.ones((self.compNo, self.compNo))

        # * check
        if res0.success is True:
            X = res0.x
            k = 0
            for i in range(self.compNo):
                for j in range(self.compNo):
                    if i != j:
                        Aij[i, j] = X[k]
                        k += 1
        else:
            return []

        #! solve nonlinear equation for alpha_ij
        fun2 = self.WilsonTemperatureIndependentParametersFunction
        # initial guess
        alpha0_ij = 0.1*np.ones(self.compNo)
        # params
        params2 = (MoVoi, Aij, T)
        res = optimize.fsolve(fun2, alpha0_ij, args=(params2,), )

        # set
        aij = np.zeros((self.compNo, self.compNo))
        k = 0
        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    aij[i, j] = res[k]
                    k += 1

        return aij

# ! NRTL

    def NRTLParameterObjectiveFunction(self, x, params):
        '''
        Non-random two-liquid model (NRTL) function

            args:
                x: parameters *** array ***
                params:
        '''
        # params
        xi_exp, ExMoGiEn_exp = params

        # parameters
        k = 0
        # temperature-dependent parameters
        taij = np.zeros((self.compNo, self.compNo))
        for i in range(self.compNo):
            for j in range(self.compNo):
                if j != i:
                    taij[i, j] = x[k]
                    k += 1

        # non-randomness parameters
        aij = np.zeros((self.compNo, self.compNo))
        for i in range(self.compNo):
            for j in range(self.compNo-i):
                if i != j+i:
                    aij[i, j+i] = x[k]
                    aij[j+i, i] = x[k]
                    k += 1

        # number of experimental data
        expDataNo = xi_exp.shape[0]

        # calculate excess molar gibbs energy
        ExMoGiEn_cal = np.zeros(expDataNo)
        for i in range(expDataNo):
            # xi
            _xi = xi_exp[i, :]

            # calculate activity coefficient
            _AcCo = self.NRTL_activity_coefficient_parameter_estimation(
                _xi, taij, aij)

            ExMoGiEn_cal[i] = self.ExcessMolarGibbsFreeEnergy(_xi, _AcCo)

        # obj function
        return ExMoGiEn_exp - ExMoGiEn_cal

    def NRTLTemperatureIndependentParametersFunction(self, x, params):
        '''
        find temperature-independent parameters (gij)

        *** gij is constant ***

        args:
            x: gij *** list ***
            params:
                1. taij: temperature dependent parameter
                2. T: temperature [K]
        '''
        # params
        taij, T = params

        # vars
        gij = np.zeros((self.compNo, self.compNo))
        y = []
        k = 0

        # set
        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    gij[i, j] = x[k]
                    k += 1

        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    _y = taij[i, j] - (gij[i, j]/(R_CONST*T))
                    y.append(_y)

        # res
        return y

    def NRTLParameterEstimator(self, params):
        '''
        parameter estimation of NRTL model for a multi-component system

        the parameters are:
            1. taij (ta[i,i]=ta[j,j]=0)
                unknownNo: component number
            2. aij: randomness parameter (a[i,j]=a[j,i])
            3. gij: interaction energy (is temperature independent)

        args:
            params:
                1. liquid mole fraction
                2. excess molar gibbs energy
                3. T: temperature
                4. boundsLimit:
                    a) taij (min, max)
                    b) aij (min, max)

        '''
        # ! check
        xi_exp, ExMoGiEn, T, boundsLimit = params

        #! least-square function for Aij
        fun1 = self.NRTLParameterObjectiveFunction

        # initial guess
        # temperature dependent parameters (taij)
        taijNo = int(self.compNo*self.compNo - self.compNo)
        # randomness parameter
        aijNo = int(taijNo/2)
        # gibbs energy of interaction between molecules
        gijNo = taijNo
        # total number
        paramsNo = taijNo + aijNo

        # params
        params1 = (xi_exp, ExMoGiEn)

        # bounds
        bU = []
        bL = []
        bounds = []
        k = 0

        # REVIEW
        A0 = np.ones(paramsNo)
        # param 1
        for i in range(taijNo):
            A0[k] = 0.1
            k += 1
            # set bounds
            # lower
            _bl = boundsLimit[0][0]
            bL.append(_bl)
            # upper
            _bu = boundsLimit[0][1]
            bU.append(_bu)

        # param 2
        for i in range(aijNo):
            A0[k] = 0.1
            k += 1
            # set bounds
            # lower
            _bl = boundsLimit[1][0]
            bL.append(_bl)
            # upper
            _bu = boundsLimit[1][1]
            bU.append(_bu)

        bounds = [bL, bU]

        res0 = optimize.least_squares(
            fun1, A0, args=(params1,), bounds=bounds)

        #! optimal NRTL parameters
        # temperature-dependent parameters
        taij = np.zeros((self.compNo, self.compNo))
        # non-randomness parameters
        aij = np.zeros((self.compNo, self.compNo))

        if res0.success is True:
            X = res0.x

            # parameters
            k = 0
            # taij
            for i in range(self.compNo):
                for j in range(self.compNo):
                    if j != i:
                        taij[i, j] = X[k]
                        k += 1
            # aij
            for i in range(self.compNo):
                for j in range(self.compNo-i):
                    if i != j+i:
                        aij[i, j+i] = X[k]
                        aij[j+i, i] = X[k]
                        k += 1
        else:
            return []

        #! solve nonlinear equation for g_ij
        fun2 = self.NRTLTemperatureIndependentParametersFunction
        # initial guess
        gij0 = 0.1*np.ones(gijNo)
        # params
        params2 = (taij, T)
        res1 = optimize.fsolve(fun2, gij0, args=(params2,), )

        # set
        gij = np.zeros((self.compNo, self.compNo))
        k = 0
        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    gij[i, j] = res1[k]
                    k += 1

        # res
        res = {
            "taij": taij,
            "aij": aij,
            "gij": gij
        }

        return res

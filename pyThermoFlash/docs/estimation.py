# PARAMETER ESTIMATION
# import libs
import numpy as np
from scipy import optimize


class Estimation:

    def __init__(self):
        pass

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

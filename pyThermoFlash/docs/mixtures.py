# MIXTURES
# import libs

# local


class Mixture:
    '''
    Mixture class to handle the mixture properties and calculations.
    '''

    def __init__(self):
        pass

    def __repr__(self):
        des = 'Mixture class to handle the mixture properties and calculations.'
        return des

        # def Txy_binary(self, pressure, guess_temperature=350, vapor_pressure_method='polynomial', model="raoult", zi_no=10):
    #     '''
    #     bubble temperature calculation

    #     args:
    #         mole_fraction: feed mole fraction (zi=xi)
    #         pressure: system pressure [Pa]
    #         guess_temperature:
    #         vapor_pressure_method: vapor-pressure calculation method (default: polynomial)
    #         model: thermodynamic model for equilibrium system (default: raoult)
    #         zi_no: number of mole fractions
    #     '''
    #     #  feed mole fraction
    #     x1 = np.linspace(0.001, 0.999, zi_no)
    #     x2 = 1 - x1
    #     #
    #     mole_fractions = np.array([x1, x2])

    #     # res
    #     _res = []

    #     for i in range(zi_no):
    #         # params
    #         params = {
    #             "zi": np.array(mole_fractions[:, i]),
    #             "P": pressure
    #         }

    #         # config
    #         config = {
    #             "Tg0": guess_temperature,
    #             "VaPeCal": vapor_pressure_method
    #         }

    #         # cal
    #         _btRes = self.bubbleTemperature(params, config)
    #         _res.append(_btRes)

    #     # Ts
    #     Ts = [item['T'] for item in _res]
    #     # xi
    #     xi = [item['xi'][0] for item in _res]
    #     # yi
    #     yi = [item['yi'][0] for item in _res]
    #     # set
    #     bubbleLine = [[xi], [Ts]]
    #     dewLine = [[yi], [Ts]]

    #     # ! plot setting
    #     # XYList = Visual.plots2DSetXYList(dataX, dataYs)
    #     # -> label
    #     # dataList = Visual.plots2DSetDataList(XYList, labelList)
    #     # Visual.plot2D(xi, Ts)
    #     # Visual.plot2D(yi, Ts)

    #     plt.plot(xi, Ts, 'r-', yi, Ts, 'b-')
    #     plt.show()

    #     # res
    #     return _res

    # def Pxy_binary(self, temperature,  vapor_pressure_method='polynomial', model="modified-raoult", activity_coefficient_model={}, zi_no=10):
    #     '''
    #     bubble temperature calculation

    #     args:
    #         mole_fraction: feed mole fraction (zi=xi)
    #         temperature: system (fixed) temperature [K]
    #         vapor_pressure_method: vapor-pressure calculation method (default: polynomial)
    #         model: thermodynamic model for equilibrium system (default: raoult)
    #         activity_coefficient_model:
    #             1. model
    #                 1.
    #             2. parameters
    #         zi_no: number of mole fractions
    #     '''
    #     #  feed mole fraction
    #     x1 = np.linspace(0.001, 0.999, zi_no)
    #     x2 = 1 - x1
    #     #
    #     mole_fractions = np.array([x1, x2])

    #     # res
    #     _res = []

    #     for i in range(zi_no):
    #         # params
    #         params = {
    #             "zi": np.array(mole_fractions[:, i]),
    #             "T": temperature
    #         }

    #         # config
    #         config = {
    #             "VaPeCal": vapor_pressure_method,
    #             "model": model,
    #             "AcCoModel": activity_coefficient_model
    #         }

    #         # ! check vle model
    #         # cal
    #         _btRes = self.bubblePressure(params, config)
    #         _res.append(_btRes)

    #     # pressure series (Ps)
    #     Ps = [item['P'] for item in _res]
    #     # xi
    #     xi = [item['xi'][0] for item in _res]
    #     # yi
    #     yi = [item['yi'][0] for item in _res]
    #     # set
    #     bubbleLine = [[xi], [Ps]]
    #     dewLine = [[yi], [Ps]]

    #     # ! plot setting
    #     # XYList = Visual.plots2DSetXYList(dataX, dataYs)
    #     # -> label
    #     # dataList = Visual.plots2DSetDataList(XYList, labelList)
    #     # Visual.plot2D(xi, Ts)
    #     # Visual.plot2D(yi, Ts)

    #     plt.plot(xi, Ps, 'r-', yi, Ps, 'b-')
    #     plt.show()

    #     # res
    #     return _res

    # def Margules_1Parameter(self, liquid_mole_fraction, vapor_mole_fraction, pressure, temperature,  vapor_pressure_method='polynomial'):
    #     '''
    #     find Margules 1-parameter based on experimental results for a binary system

    #     args:
    #         liquid_mole_fraction: liquid mole fraction as [x1,x2]
    #         vapor_mole_fraction: vapor mole fraction as [y1,y2]
    #         pressure: system pressure as P

    #     return:
    #         A12: Margules 1-parameter
    #     '''
    #     try:
    #         # vapor pressure at inlet temperature
    #         VaPr = self.vaporPressureMixture(
    #             temperature, vapor_pressure_method)

    #         # calculate excess molar gibbs energy
    #         #! check
    #         if len(liquid_mole_fraction) == 1 and len(vapor_mole_fraction) == 1:
    #             # activity coefficient from experimental data
    #             AcCoExp = self.activityCoefficientUsingModifiedRaoult(
    #                 liquid_mole_fraction[0], vapor_mole_fraction[0], pressure[0], VaPr)

    #             # one algebraic equation (A(12)*x1*x2)
    #             ExMoGiEn = self.ExcessMolarGibbsFreeEnergy(
    #                 liquid_mole_fraction[0], AcCoExp)

    #             # A12 parameter
    #             A12 = ExMoGiEn / \
    #                 (liquid_mole_fraction[0][0]*liquid_mole_fraction[0][1])

    #             return A12
    #         else:
    #             # solve a system of non-linear equation
    #             pass

    #     except Exception as e:
    #         pass

    # def Margules_parameter_estimation(self, csv_file, mode="2-parameter", vapor_pressure_method='polynomial', plot_result=True):
    #     '''
    #     estimate Margules parameters for a *** binary system ***

    #     args:
    #         csv_file: csv file path with data format as:
    #             1. no
    #             2. T: fixed temperature [K]
    #             3. P: measured pressure [Pa]
    #             4. x1: liquid mole fraction component 1
    #             5. y1: vapor mole fraction component 1

    #     '''
    #     np_data, df_data, rowNo, colNo, colsName = LoaddataClass.load_csv_to_df(
    #         csv_file)

    #     # parameters
    #     parameterNo = 0

    #     # ! check
    #     if mode == '2-parameter':
    #         # set
    #         parameterNo = 2
    #     else:
    #         # set
    #         parameterNo = 1

    #     # calculate activity coefficient using modified-raoult's law
    #     ExMoGiEn = np.zeros(rowNo)

    #     # interpret Pxy data (binary system)
    #     AcCo, xi = LoaddataClass.Pxy_BinarySystemInterpretData(
    #         self.pool, np_data, rowNo, vapor_pressure_method)

    #     for i in range(rowNo):
    #         # check for x1=0, activity coefficient is not defined
    #         if i > 0 and i < rowNo-1:
    #             # calculate excess molar gibbs energy
    #             ExMoGiEn[i] = self.ExcessMolarGibbsFreeEnergy(
    #                 xi[i, :], AcCo[i, :])

    #     #! call optimizer fun
    #     # params
    #     params = (xi[1:-1, :], ExMoGiEn[1:-1], parameterNo)

    #     res = self.margulesParameterEstimator(params)

    #     if res.success is True:
    #         Aij = res.x

    #         # use for the model
    #         ExMoGiEn_model = np.zeros(rowNo)
    #         AcCo_model = np.zeros((rowNo, 2))

    #         for i in range(rowNo):
    #             if i > 0 and i < rowNo-1:
    #                 # activity coefficient using the model
    #                 AcCo_model[i, :] = self.Margules_activity_coefficient(
    #                     xi[i, :], Aij)

    #                 ExMoGiEn_model[i] = self.ExcessMolarGibbsFreeEnergy(
    #                     xi[i, :], AcCo_model[i, :])

    #     # plot
    #     if plot_result is True:
    #         plt.plot(xi[1:-1, 0], ExMoGiEn[1:-1], 'o')
    #         plt.plot(xi[1:-1, 0], ExMoGiEn_model[1:-1], 'b-')
    #         plt.show()

    #     # res
    #     return Aij

    # def Wilson_parameter_estimation(self, csv_file, vapor_pressure_method='polynomial', plot_result=True, aij_bounds=[0, 10]):
    #     '''
    #     estimate wilson parameters for a *** binary system ***

    #     args:
    #         csv_file: csv file path with data format as:
    #             1. no
    #             2. T: fixed temperature [K]
    #             3. P: measured pressure [Pa]
    #             4. x1: liquid mole fraction component 1
    #             5. y1: vapor mole fraction component 1
    #         vapor_pressure_method: vapor pressure method
    #         plot_result: to display results (default: true)
    #         aij_bounds: alpha[i,j] lower and upper bounds

    #     return:
    #         alpha: composition-independent parameters
    #         that describe how the interactions between the unlike components differ from the like components

    #     '''
    #     np_data, df_data, rowNo, colNo, colsName = LoaddataClass.load_csv_to_df(
    #         csv_file)

    #     # calculate activity coefficient using modified-raoult's law
    #     ExMoGiEn = np.zeros(rowNo)

    #     # interpret Pxy data (binary system)
    #     AcCo, xi, T = LoaddataClass.Pxy_BinarySystemInterpretData(
    #         self.pool, np_data, rowNo, vapor_pressure_method)

    #     for i in range(rowNo):
    #         # check for x1=0, activity coefficient is not defined
    #         if i > 0 and i < rowNo-1:
    #             # calculate excess molar gibbs energy
    #             ExMoGiEn[i] = self.ExcessMolarGibbsFreeEnergy(
    #                 xi[i, :], AcCo[i, :])

    #     #! call optimizer fun
    #     # params
    #     params = (xi[1:-1, :], ExMoGiEn[1:-1], T)

    #     aij = self.WilsonParameterEstimator(params)

    #     # REVIEW
    #     # use for the model
    #     ExMoGiEn_model = np.zeros(rowNo)
    #     AcCo_model = np.zeros((rowNo, 2))

    #     for i in range(rowNo):
    #         if i > 0 and i < rowNo-1:
    #             # activity coefficient using the model
    #             AcCo_model[i, :] = self.Wilson_activity_coefficient(
    #                 xi[i, :], T, aij)

    #             ExMoGiEn_model[i] = self.ExcessMolarGibbsFreeEnergy(
    #                 xi[i, :], AcCo_model[i, :])

    #     # plot
    #     if plot_result is True:
    #         plt.plot(xi[1:-1, 0], ExMoGiEn[1:-1], 'o')
    #         plt.plot(xi[1:-1, 0], ExMoGiEn_model[1:-1], 'b-')
    #         plt.show()

    #     # res
    #     return aij

    # def NRTL_parameter_estimation(self, csv_file, vapor_pressure_method='polynomial', plot_result=True, bounds=[[-5, 5], [-5, 5]]):
    #     '''
    #     estimate NRTL parameters for a *** multi-component system ***

    #     args:
    #         csv_file: csv file path with data format as:
    #             1. no
    #             2. T: fixed temperature [K]
    #             3. P: measured pressure [Pa]
    #             4. x1: liquid mole fraction component 1
    #             5. y1: vapor mole fraction component 1
    #         vapor_pressure_method: vapor pressure method
    #         plot_result: to display results (default: true)
    #         bounds: lower and upper bounds for
    #             a) taij: temperature dependent parameter
    #             b) aij: randomness parameter
    #             c) gij: interaction energy parameter

    #     return:
    #         taij: temperature dependent parameter
    #         aij: randomness parameter
    #         gij: interaction energy parameter

    #     '''
    #     np_data, df_data, rowNo, colNo, colsName = LoaddataClass.load_csv_to_df(
    #         csv_file)

    #     # calculate activity coefficient using modified-raoult's law
    #     ExMoGiEn = np.zeros(rowNo)

    #     # interpret Pxy data (binary system)
    #     AcCo, xi, T = LoaddataClass.Pxy_BinarySystemInterpretData(
    #         self.pool, np_data, rowNo, vapor_pressure_method)

    #     for i in range(rowNo):
    #         # check for x1=0, activity coefficient is not defined
    #         if i > 0 and i < rowNo-1:
    #             # calculate excess molar gibbs energy
    #             ExMoGiEn[i] = self.ExcessMolarGibbsFreeEnergy(
    #                 xi[i, :], AcCo[i, :])

    #     #! call optimizer fun
    #     # params
    #     params = (xi[1:-1, :], ExMoGiEn[1:-1], T, bounds)
    #     res = self.NRTLParameterEstimator(params)

    #     # set
    #     taij = res['taij']
    #     aij = res['aij']
    #     gij = res['gij']

    #     # REVIEW
    #     # use for the model
    #     ExMoGiEn_model = np.zeros(rowNo)
    #     AcCo_model = np.zeros((rowNo, 2))

    #     for i in range(rowNo):
    #         if i > 0 and i < rowNo-1:
    #             # activity coefficient using the model
    #             AcCo_model[i, :] = self.NRTL_activity_coefficient(
    #                 xi[i, :], T, aij, gij)

    #             ExMoGiEn_model[i] = self.ExcessMolarGibbsFreeEnergy(
    #                 xi[i, :], AcCo_model[i, :])

    #     # plot
    #     if plot_result is True:
    #         plt.plot(xi[1:-1, 0], ExMoGiEn[1:-1], 'o')
    #         plt.plot(xi[1:-1, 0], ExMoGiEn_model[1:-1], 'b-')
    #         plt.show()

    #     # res
    #     return res

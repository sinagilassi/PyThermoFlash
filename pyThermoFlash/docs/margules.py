# MARGULES MODELS
# ----------------

# packages/modules
from math import pow, sqrt, exp
import numpy as np
# local


class Margules:
    def __init__(self, components):
        self.components = components
        # set
        self.compNo = len(components)

    def Margules_activity_coefficient(self, xi, Aij):
        '''
        the Gibbs excess function is defined as a one-parameter Margules equation (binary system)

        args:
            xi: mole fraction
            Aij: Margules constant (which is taken from experimental results)
        '''
        # check
        if self.compNo > 2:
            return [0, 0]

        # activity coefficient
        AcCo = []

        # ! check method
        if Aij.size == 1:
            # 1-parameter
            for i in range(self.compNo):
                AcCo.append(exp(Aij[0]*pow(1-xi[i], 2)))
        elif Aij.size == 2:
            # 2-parameter
            # set
            A12 = Aij[0]
            A21 = Aij[1]
            for i in range(self.compNo):
                if i == 0:
                    _AcCo = exp(pow(1-xi[i], 2) *
                                (A12 + 2*(A21 - A12)*xi[i]))
                elif i == 1:
                    _AcCo = exp(pow(1-xi[i], 2) *
                                (A21 + 2*(A12 - A21)*xi[i]))

                AcCo.append(_AcCo)

        # res
        return AcCo

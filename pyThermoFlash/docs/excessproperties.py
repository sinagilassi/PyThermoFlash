# EXCESS PROPERTIES
# ------------------

# packages/modules
from math import pow, sqrt, exp, log
import numpy as np
# local
from PyCTPM.core import R_CONST


class ExcessProperties:
    def __init__(self, components):
        self.components = components
        # set
        self.compNo = len(components)

    def VanLaar_am_bm(self, xi, ai, bi):
        '''
        van der Waals one-fluid mixing rules

        return:
            am
            bm
        '''
        # a
        am = []
        bm = []
        for i in range(self.compNo):
            _xi = xi[i]
            _ai = ai[i]
            _bi = bi[i]
            for j in range(self.compNo):
                # set parameter
                if i == j:
                    _aij = ai[j]
                    _bij = bi[j]
                else:
                    _aij = sqrt(_ai*ai[j])
                    _bij = (_bi + bi[j])/2
                # set
                _ami = _xi*xi[j]*_aij
                am.append(_ami)
                _bmi = _xi*xi[j]*_bij
                bm.append(_bmi)

        # sum
        am = np.array(am).sum()
        bm = np.array(bm).sum()

        # res
        return am, bm

    def VanLaar_L(self, ai, bi, T):
        '''
        calculate L[i,j]
        '''
        # L
        Lij = np.zeros((self.compNo, self.compNo))

        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    _c0 = sqrt(ai[i])/bi[i]
                    _c1 = sqrt(ai[j])/bi[j]
                    _c2 = pow(_c0-_c1, 2)
                    Lij[i, j] = (bi[i]/(R_CONST*T))*_c2

        # res
        return Lij

    def VanLaar_activity_coefficient(self, xi, ai, bi, T):
        '''
        calculate activity coefficient

        args:
            xi: mole fraction 
            ai: van der Waals a constant
            bi: van der Waals b constant 
        '''
        # gamma
        AcCo = np.zeros(self.compNo)

        # am/bm
        # am, bm = self.VanLaar_am_bm(ai, bi)

        # L
        Lij = self.VanLaar_L(ai, bi, T)

        for i in range(self.compNo):
            for j in range(self.compNo):
                if i != j:
                    _res0 = Lij[i, j] * \
                        pow(1+((Lij[i, j]*xi[i])/(Lij[j, i]*xi[j])), -2)
                    _res1 = exp(_res0)
                    # res
                    AcCo[i] = _res1

        # res
        return AcCo

    def ExcessMolarGibbsFreeEnergy(self, xi, AcCoi):
        '''
        calculate excess molar Gibbs free energy from activity coefficient of each component in a system

        args:
            xi: liquid mole fraction
            AcCoi: activity coefficient
        '''
        # G(E)/RT
        c0 = np.zeros(self.compNo)

        for i in range(self.compNo):
            c0[i] = xi[i]*log(AcCoi[i])

        # excess molar gibbs energy
        c1 = np.sum(c0)

        # res
        return c1

# DISPLAY RESULTS
# ----------------

# packages/modules
import numpy as np


class Display:
    def __init__(self):
        pass

    def colDisplay(self, data):
        colCenterLen = 50
        # set
        dash = '-' * (colCenterLen*2)

        # col length set
        colLenManualSet = round(colCenterLen/2)

        for i in range(len(data)):
            if i == 0:
                print(dash)
                print('{:^{colCenterLen}s}{:^{colLenManualSet}s}{:^{colLenManualSet}s}'.format(
                    data[i][0], data[i][1], data[i][2], colCenterLen=colCenterLen, colLenManualSet=colLenManualSet))
                print(dash)
            else:
                print('{:<{colCenterLen}s}{:<{colLenManualSet}}{:>0}'.format(
                    data[i][0].upper(), data[i][1], data[i][2], colCenterLen=colCenterLen, colLenManualSet=colLenManualSet))

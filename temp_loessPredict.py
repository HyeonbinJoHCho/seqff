import numpy as np
import matplotlib.pyplot as plt
from skmisc.loess import loess

def loess_debugging(loessVar):
    y = loessVar.predict(np.array([0.5, 0.501, 0.505, 0.510]))
    # y = loess.loess_prediction([0.5, 0.7], loessVar)    # Extrapolation not allowed
    # print(y)

    print(y.values)
    print(loessVar.model)

    # loessVar.model.span = 0.10                 #change loessVar span
    # loessVar.fit()
    # print(loessVar.model)

    xTemp = np.arange(0.317, 0.651, 0.001)
    yTemp = loessVar.predict(xTemp)

    plt.hist2d(xTemp, yTemp.values, bins=len(xTemp))
    plt.show()



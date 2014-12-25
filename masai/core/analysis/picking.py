import numpy as np
import nmrglue as ng


def pick(source, thres=None, index=False):
    if thres is None:
        thres = 80.

    pthres = np.abs(source.data.values).max() * thres / 100.

    ar = ng.analysis.peakpick.pick(np.abs(source.data.values), pthres, table=False)[0]

    if not index:
        x = source.data.columns
        y = source.data.index
        ars = []
        for i, item in enumerate(ar):
            ar1 = x[item[1]]
            ar0 = y[item[0]]
            ars.append([ar0, ar1])
        return ars
    else:
        return ar



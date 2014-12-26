import numpy as np
from nmrglue.analysis.peakpick import pick


#
# picking
#
def picking(source, thres=None, index=False):
    """

    Parameters
    ----------
    source: source instance
    thres: float, optional
        percent of threshold
    index:

    Returns
    -------
    """
    if thres is None:
        thres = 80.

    pthres = np.abs(source.data.values).max() * thres / 100.

    ar = pick(np.abs(source.data.values), pthres, table=False)[0]

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



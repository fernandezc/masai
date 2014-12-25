__author__ = 'christian'

# ===============================================================================
# process_stmas
#===============================================================================
def process_stmas(source, q=2, **kargs):
    """
    STMAS processing

    """
    kargs['shearing'] = (1, q)
    print('process_stmas... ')
    return _processstmq(source, **kargs)


#===============================================================================
# process_mqmas
#===============================================================================
def process_mqmas(source, p=3, **kargs):
    """
    MQMAS processing

    """
    kargs['shearing'] = (p, 0)
    print('process_mqmas... ')

    return _processstmq(source, **kargs)


# =======================================================================================================================
# linear prediction
#=======================================================================================================================
def lp(source, axis=-1, **kargs):
    if kargs.get('debug', False):
        print("linear prediction")
    if axis == 0:
        k = 'lpbin1'
        data = source.data.values.T.copy()
    else:
        k = 'lpbin'
        data = source.data.values.copy()
    pred = source.fetchPar('LPBIN', axis=axis)
    pred = kargs.get(k, pred)
    data = linear_predict(data, pred)
    if data.shape[-1] % 2 != 0:
        data = data[:, :-1]  # make it even
    source.storePar('TD', data.shape[-1], axis=axis)
    source.storePar('TDeff', data.shape[-1], axis=axis)
    if axis == 0:
        source.data = data.T
    else:
        source.data = data

#####################
# private functions #
#####################

#=======================================================================================================================
# processstmq
#=======================================================================================================================
def _processstmq(source, **kargs):
    """
    MQMAS processing

    """
    basemode = kargs.get('basemode', 'basf')
    print(basemode)
    xf2(source, **kargs)
    if basemode is not None:
        autobase(source, basetype=basemode, **kargs)
    xf1(source, **kargs)
    if kargs.get('lastft', True) and (basemode is not None):
        autobase(source, basetype=basemode, axis=0, **kargs)


#===============================================================================
# makestates
#===============================================================================
def _makestates(source, tppi=False):
    """
    Transforms according to the States or states-tppi procedure.

    """
    data = source.data.values

    ar = data[::2]  # transposed real part
    ai = data[1::2]  # transposed imag part
    data = ar.real + ai.real * 1.0j
    #data = ar.imag + ai.imag * 1.0j

    if tppi:
         data[1::2] = -data[1::2]

    source.storePar('TD', data.shape[0], axis=0)

    source.data = data

#===============================================================================
# _apodization
#===============================================================================
def _apodization(source, func, lb, top, axis=-1):
    """
    Apply a given apodization to the source data

    """
    model = ApodModels(source=source)

    _func = getattr(model, func)

    if axis == 0 and source.is_2D:
        source.data[:] = (source.data.values.T * _func(lb, top, 0)).T
    else:
        source.data[:] = (source.data.values * _func(lb, top, -1))


####################################################################################
class ApodModels(HasTraits):
    """
    A class to handle apodisation models

    """

    source = Instance(NMRSource)

    #===============================================================================
    # em_func
    #===============================================================================
    def em_func(self, lb, top, axis=-1):
        """
        Calculate an exponential apodization fonction

        Parameters
        ----------
        lb: float
            apodization value in spectral units
        top: float
            k shift ratio
        axis : optional, default is -1

        Returns
        -------
        f : ndarray
            the apodisation array

        """

        data = self.source.data
        lb = lb / self.source.getMultiplicator(axis)

        if axis == 0 and self.source.is_2D:
            x = data.index.values
        else:
            x = data.columns.values

        return np.exp(-np.pi * np.abs(x - top) * lb)

    #===============================================================================
    # k_shifted_em_func
    #===============================================================================
    def k_shifted_em_func(self, lb, top, axis=-1):
        """
        Calculate a shifted exponential apodization fonction

        Parameters
        ----------
        lb: float
            apodization value in spectral units
        top: float
            k shift ratio
        axis : optional, default is -1

        Returns
        -------
        f : ndarray
            the apodisation array

        """
        if not self.source.is_2D:
            raise Exception('only work for two dimensional data')
        data = self.source.data
        lb = lb / self.source.getMultiplicator(axis)

        if axis == 0:
            data = data.T

        x = data.columns.values
        y = data.index.values

        ratio = abs(top)
        # ratio must be converted in points ratio
        r = (y[1]-y[0])/(x[1]-x[0])
        ratio = ratio * r

        f = np.fromfunction(lambda y, x, lb, ratio: np.exp(-np.pi * np.abs(x - y*ratio) * lb),
                            data.shape, dtype=float, lb=lb, ratio=ratio)
        # If we are in states mode (general case two rows correspond
        # to the same time: we need to correct f for this
        # TODO: check this from parameters

        if axis == -1:
            if self.source.fnmode.startswith('STATES'):
                f[1::2] = f[::2]

        return f


# TODO: other apodization function
# from nmrglue :

# gm(data[, g1, g2, g3, inv, rev]) 	Lorentz-to-Gauss apodization
#gmb(data[, a, b, inv, rev]) 	Modified gaussian apodization
#jmod(data[, e, off, end, inv, rev]) 	Exponentially damped J-modulation apodization
#sp(data[, off, end, pow, inv, rev]) 	Shifted sine-bell apodization
#sine(data[, off, end, pow, inv, rev]) 	Shifted sine-bell apodization
#tm(data[, t1, t2, inv, rev]) 	Trapezoid Apodization
#tri(data[, loc, lHi, rHi, inv, rev]) 	Triangle apodization.

# ===============================================================================
# MAIN
#===============================================================================
if __name__ == "__main__":
    # import doctest
    #doctest.testmod()

    import os
    import numpy as np

    from lcslib.api import *

    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))

    # 1D

    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))

    fd = NMRSource(path=path)

    # undersampling
    undersampling(fd, 1)

    # linebroadening
    store(fd, 'LB', 30.)
    em(fd)

    # fourier transform
    store(fd, 'SI', 4096)
    ft(fd)

    # clipping
    clip(fd, -50, 70)

    # automatic phasing
    store(fd, 'PHC0', 0.)
    store(fd, 'PHC1', 0.)
    apk(fd, optmethod='simplex', apkmode='negmin+entropy', gamma=1., use_phc1=False, prints=True, debug=True)

    # plot
    p = Plotting(fd)
    p.figsize = (6., 3.)
    ax = p.plot(showz=True, xlim=[-45, 65], xreverse=True, color='blue')

    # automatic baseline correction
    autobase(fd, window=.025)
    autobase(fd, basetype='basf', nbzone=32, mult=5, multpolynom=True, porder=6, apply=True)
    ax = p.plot(showz=True, xlim=[-45, 65], xreverse=True, color='blue')

    ax.plot(fd.data.columns, fd.real.T, ':')
    ax.axhline(0, color='c', linestyle='--')
    p.show()



    # 2D
    user = 'aurelie'
    name = 'RbNO3'
    expno = 15
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))

    fd = NMRSource(path=path)

    undersampling(fd, 2)
    clip(fd, 0, 3.68, axis=0)

    store(fd, 'LB', 5.)
    p = 1
    q = 2
    store(fd, 'coherence_order', p)  # coherence for stmas (both p and q required with k_shifted method
    store(fd, 'satellite_order', q)  # satellite order for stmas

    store(fd, 'PHC1', 0)
    store(fd, 'PHC1', 0, axis=0)
    store(fd, 'PHC0', 0)
    store(fd, 'PHC0', 0, axis=0)

    process_stmas(fd, method='k_shifted',
                  apk='first', apk1='col', thres=50,
                  use_phc1=True, bound_phc1=10.0,
                  apkmode='negmin+entropy', entropyd=1, gamma=.1,
                  prints=True,
                  optmode='hopping',  # 'simplex',  # 'hopping'
                  basemode='linear',
                  mult=5,
                  lp1=True,
                  lpbin1=128,
                  debug=True)

    p = Plotting(fd)
    p.colormap = 'jet'
    p.xreverse = True
    p.proj = 'sky'
    p.start = -1
    p.contour()

    print("*** kill plot to exit ***")
    p.show()



from __future__ import print_function, division

#
# imports
#
import numpy as np
from nmrglue.process.proc_base import ps
from warnings import warn

from masai.core.fit.optimize import optimize, external
from masai.core.process.baseline import basecorr
from masai.core.analysis.picking import picking


#
# apk
#
def apk(source, axis=-1, **kargs):
    """
    Automatic Phasing

    Parameters
    ----------
    axis: optional, default -1
    apkmode: str, optional
        first, max, col, ... useful for 2D data

    """
    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if axis != 0:
        _apk = kargs.get('auto', 'first')
    else:
        _apk = kargs.get('auto', 'max')

    kargs['auto'] = _apk

    pk(source, axis=axis, **kargs)


#
# pk
#
def pk(source, axis=-1, **kargs):
    """
    Automatic or manual phasing

    Parameters
    ----------
    use_phc1: bool, optional, default=False
        also optimize first order phase if True
    mode: String, optional, default='negmin'
        method for automatic phase detection: 'negmin', 'entropy', ...
    entropyd: int, optional, default=2
        order of derivation for method entropy
    gamma: float, optional
        relative weigth for method entropy error with respect to negmin
    axis: optional, default=-1

    """

    if axis == -1:
        par = source.par
    else:
        par = source.par2

    if not par.isfreq:
        warn('This source is not already transformed: are you sure this is what you want?')

    # we will work on the ndarray, not DataFrame
    data = source.data.values

    if axis == 0:
        # transpose temporarily the data for indirect dimension ft
        data = data.T

    # get the initial phase setting
    phc0 = par.PHC0
    phc1 = par.PHC1

    if not kargs.get('auto'):
        # not automatic leaves here...
        data = ps(data, phc0, phc1)

    else:
        # APK
        if kargs.get('auto') == 'col':
            ar = picking(source, kargs.get('thres', 50), index=True)
            kargs['cols'] = zip(*ar)[0]
        data, phc0, phc1 = autophase(data, phc0, phc1, **kargs)

        # store the new phases
        par.PHC0 = phc0
        par.PHC1 = phc1

    # un-transpose the data if needed
    if axis == 0:
        data = data.T

    source.data = data

    # apply linear baseline correction
    # autobase(source, axis=axis)


#
# autophase
#
def autophase(data, phc0=0.0, phc1=0.0, **kargs):
    """
    Automatic phasing of 1D or 2D spectra

    Parameters
    ----------
    data: ndarray
        Complex data to be phased
    phc0: float, optional
        Zero-order phase
    phc1: float, optional
        First-order phase
    apk: any, optional, default is 'first'
        If 'first', only the first raw is taken into account for phasing
        If 'proj', phase is calculated on the projection
    apkmode: string, optional
        Method of optimization among 'negmin' or 'entropy'. Default is 'negmin+entropy'
    entropyd: int, optional, default is 1
        Order of the derivative to use for the calculation of the entropy
    gamma: float, optional, default is .1
        Weight of the 'entropy' optimization relative to the 'negmin' optimization
    print: bool, optional, default is False
      If True some results information are printed

    Returns
    -------
    data: ndarray
        phased data
    phc0: float
        New phc0 phase
    phc1: float
        New phc1 phase

    """
    dat = data.copy()

    # get keyword parameters
    # ------------------------
    debug = kargs.get('debug', False)
    printflag = kargs.get('prints', False)

    if debug:
        print('performing autophase')

    apk = kargs.get('auto', 'first')

    if apk in ['first', 'standard'] or data.shape[0] == 1:
        # select the first row
        dat = dat[0]
        if data.shape[0] == 1:
            apk = 'standard'

    if apk == 'max' and data.ndim > 1:
        i, j = np.unravel_index(data.real.argmax(), data.shape)
        dat = data[i]

    if apk == 'col' and data.ndim > 1:
        l = []
        cols = kargs['cols']
        if printflag:
            print('columns selected:', cols)
        for col in cols:
            i = int(col)
            l.append(data[i])
        dat = np.vstack(l)

    # apkmode
    if debug:
        print('\napk: ', apk)

    apkmode = kargs.get('apkmode', 'negmin+entropy')
    entropyd = kargs.get('entropyd', 1)
    gamma = kargs.get('gamma', .1)
    optmode = kargs.get('optmode', 'SIMPLEX')

    # baseline mode
    basemode = kargs.get('basemode', 'linear')
    if basemode == 'basf':
        print('warning, basemode=basf should be avoided. This can take a very long time')

    # prepare the phase to optimized with bound
    bp0 = kargs.get('bound_phc0', 360.)
    bp1 = kargs.get('bound_phc1', 10.)

    p0s = [phc0, -bp0, bp0]
    bounds = [p0s[1:], ]
    p = [p0s, ]

    use_phc1 = kargs.get('use_phc1', False)
    if use_phc1:
        p1s = [phc1, -bp1, bp1]
        p.append(p1s)
        bounds.append(p1s[1:])

    # optim
    global ni, nas, spe, err
    ni = 0
    err = 0

    # internal error function ------------------------------------------------------------------------------------------
    def _phase_error(p, bounds, s):

        global nas, spe, ni, err

        # restore the external parameter
        pe = []
        for i, x in enumerate(bounds):
            if not isinstance(x, list):
                pe.append(p[i])
            else:
                # restore external parameters in case of bounding
                pe.append(external(p[i], x[0], x[1]))

        p0 = pe[0]

        sc = s.copy()
        p1 = 0.0
        if use_phc1:
            p1 = pe[1]
            sc = ps(sc, p0, p1)
        else:
            sc = ps(sc, p0)

        # baseline correction
        scp = sc - basecorr(sc, basetype=basemode, **kargs)

        # Negative area minimization
        nas = 1.
        nam = 1.
        if "negmin" in apkmode:
            fm = scp[scp.real <= 0].real
            fm = fm - np.min(fm)
            nas = np.sum(fm ** 2)

            # normalisation
            fp = scp.real
            fp = fp - np.min(fm)
            nas = nas * 100. / (np.sum(fp ** 2) + 1.e-30)

        # entropy minimization
        spe = 1.
        if 'entropy' in apkmode:
            h = np.diff(scp.real, entropyd)
            h = np.abs(h)
            h = abs(h / np.sum(np.abs(h))) + 1.e-30
            spe = -np.sum(h * np.log(h))
        ni += 1
        err = nas + spe * gamma

        return err

    # end _phase_error function -------------------------------------------------------------------------------------

    # callback function--------------------------------------------------------

    def print_fun(x, f=None, accepted=None):
        """
        callback print function
        """
        global nas, spe, ni, err
        if debug and f is not None:  # hopping
            print("i: %d >>> NA: %.3g, S: %.3g err: %.3g (accepted %d)" % (
                ni, nas, spe, f, int(accepted)))
        if debug and f is None:  # simplex
            print("i: %d >>> NA: %.3g, S: %.3g err: %.3g" % (
                ni, nas, spe, err))

    # end callback function ---------------------------------------------------

    # convergence is not insured depending on the starting values
    phase, err1 = optimize(_phase_error, p,
                           args=(bounds, dat), method=optmode,
                           callback=print_fun)
    nasave = nas
    ssave = spe
    p = phase[:]

    p0s = [(phase[0] - 180.) % 360.0, -bp0, bp0]
    phase1 = [p0s, ]
    if len(phase) > 1:
        p1s = [phase[1], -bp1, bp1]
        phase1.append(p1s)

    if optmode.upper() != 'HOPPING':
        ni = 0
        phase2, err2 = optimize(_phase_error, phase1,
                                args=(bounds, dat), method=optmode,
                                callback=print_fun)

        # select the best
        p = phase2[:]
        if err2 > err1:
            p = phase[:]
            nas = nasave
            spe = ssave

    # extract results

    if printflag:
        print('Phase:')
    if use_phc1:
        phc0, phc1 = p
        if printflag:
            print('phc0', phc0)
        if printflag:
            print('phc1', phc1)
    else:
        phc0 = p[0]
        phc1 = 0.0
        if printflag:
            print('phc0', phc0)
        if printflag:
            print('phc1', phc1, 'not optimized')

    if printflag:
        print('Negative area (NA):%.3g Entropy (S) * gamma:%.3g' % (nas, spe * gamma))

    # apply to the original data and return
    data = ps(data, phc0, phc1)
    return data, phc0, phc1


if __name__ == '__main__':
    optimize()

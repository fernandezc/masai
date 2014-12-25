# =======================================================================================================================
# autophase
# =======================================================================================================================
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
    #------------------------
    debug = kargs.get('debug', False)
    printflag = kargs.get('prints', False)

    if debug: print('performing autophase')

    apk = kargs.get('apk','first')

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
    if debug: print('\napk: ', apk)
    apkmode = kargs.get('apkmode','negmin + entropy')
    entropyd = kargs.get('entropyd',1)
    gamma = kargs.get('gamma', .1)
    optmode = kargs.get('optmode', 'SIMPLEX')

    # baseline mode
    basemode = kargs.get('basemode','linear')
    if basemode == 'basf':
        print ('warning, basemode=basf should be avoided. This can take a very long time')

    # prepare the phase to optimized with bound
    bp0 = kargs.get('bound_phc0',360.)
    bp1 = kargs.get('bound_phc1', 10.)

    p0s = [phc0, -bp0, bp0]
    bounds = [p0s[1:],]
    p =[p0s,]

    use_phc1 = kargs.get('use_phc1', False)
    if use_phc1:
        p1s = [phc1, -bp1, bp1]
        p.append(p1s)
        bounds.append(p1s[1:])

    # optimize
    global ni, NA, S, err
    ni = 0

    # internal error function ------------------------------------------------------------------------------------------
    def _phaseError(p, bounds, s):

        global NA, S, ni, err

        #restore the external parameter
        pe = []
        for i, x in enumerate(bounds):
            if not isinstance(x, list):
                pe.append(p[i])
            else:
                # restore external parameters in case of bounding
                pe.append(external(p[i], x[0], x[1]))

        p0 = pe[0]

        sc = s.copy()
        p1=0.0
        if use_phc1:
            p1 = pe[1]
            sc = ng.process.proc_base.ps(sc, p0, p1)
        else:
            sc = ng.process.proc_base.ps(sc, p0)

        # baseline correction
        scp = sc - basecorr(sc, basetype=basemode, **kargs )

        # Negative area minimization
        NA = 1.
        NAm = 1.
        if "negmin" in apkmode:
            fm = scp[scp.real<=0].real
            fm = fm - np.min(fm)
            NA = np.sum(fm**2)

            # normalisation
            fp = scp.real
            fp = fp - np.min(fm)
            NA = NA *100. / (np.sum(fp**2)+1.e-30)

        # entropy minimization
        S = 1.
        if 'entropy' in apkmode:
            h = np.diff(scp.real,entropyd)
            h = np.abs(h)
            h = abs(h / np.sum(np.abs(h))) + 1.e-30
            S = -np.sum(h*np.log(h))
        ni += 1
        err = NA + S*gamma

        return err

    # end _phaseerror function -------------------------------------------------------------------------------------

    # callback function--------------------------------------------------------

    def print_fun(x, f=None, accepted=None):
        """
        callback print function
        """
        global NA, S, ni, err
        if debug and f is not None:  # hopping
            print("i: %d >>> NA: %.3g, S: %.3g err: %.3g (accepted %d)" % (
                ni, NA, S, f, int(accepted)))
        if debug and f is None:  # simplex
            print("i: %d >>> NA: %.3g, S: %.3g err: %.3g" % (
                ni, NA, S, err))

    # end callback function ---------------------------------------------------

    # convergence is not insured depending on the starting values
    phase, err1 = optimize(_phaseError, p,
                           args=(bounds, dat), method=optmode,
                          callback=print_fun)
    NAsave = NA
    Ssave = S
    p = phase[:]

    p0s =[(phase[0]-180.)%360.0, -bp0, bp0]
    phase1 = [p0s,]
    if len(phase)>1:
        p1s =[phase[1], -bp1, bp1]
        phase1.append(p1s)

    if optmode.upper() != 'HOPPING':
        ni = 0
        phase2, err2 = optimize(_phaseError, phase1,
                                args=(bounds, dat), method=optmode,
                                callback=print_fun)

        # select the best
        p = phase2[:]
        if err2 > err1:
            p = phase[:]
            NA = NAsave
            S = Ssave

    # extract results

    if printflag: print('Phase:')
    if use_phc1:
        phc0, phc1 = p
        if printflag: print('phc0', phc0)
        if printflag: print('phc1', phc1)
    else:
        phc0 = p[0]
        if printflag: print('phc0', phc0 )
        if printflag: print('phc1', phc1, 'not optimized')

    if printflag: print('Negative area (NA):%.3g Entropy (S) * gamma:%.3g'%(NA, S*gamma))

    #apply to the original data and return
    data = ng.process.proc_base.ps(data, phc0, phc1)
    return data, phc0, phc1

#===============================================================================
# pk
# ===============================================================================
def pk(source, axis=-1, **kargs):
    """
    Automatic or manual phasing

    Parameters
    ----------
    apk: bool, optional
        automatic phasing if True
        * for 2D spectra:
            * if True, the phasing is automatic on all rows
            * if first, only the first row is phased and the phase applied to the whole 2D spectrum
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

    # automatic phasing?
    apk = False
    apk = kargs.get('apk', None)

    # we will work on the ndarray, not DataFrame
    data = source.data.values

    # get the initial phase setting
    phc0 = source.fetchPar('PHC0', axis=axis)
    phc1 = source.fetchPar('PHC1', axis=axis)

    if apk is None:
        # not automatic leaves here...
        source.history.append('pk %s' % axis)
        if axis == 0 and source.is_2D:
            data = ng.process.proc_base.ps(data.T, phc0, phc1).T
        else:
            data = ng.process.proc_base.ps(data, phc0, phc1)

    else:
        # APK -----
        source.history.append('apk %s' % axis)

        if axis == -1 or axis == 1:
            data, phc0, phc1 = autophase(data, phc0, phc1, **kargs)
        else:
            if apk == 'col':
                ar = pick(source, kargs.get('thres', 50), index=True)
                kargs['cols'] = zip(*ar)[0]
            data, phc0, phc1 = autophase(data.T, phc0, phc1, **kargs)
            data = data.T

        # store the new phases
        source.storePar('PHC0', phc0, axis=axis)
        source.storePar('PHC1', phc1, axis=axis)

    source.data[:] = data

    # apply linear baseline correction
    autobase(source)


#===============================================================================
# apk
#===============================================================================
def apk(source, axis=-1, **kargs):
    """
    Automatic Phasing

    Parameters
    ----------
    axis: optional, default -1
    kargs: optional arguments, see 'pk'

    """
    if axis != 0:
        apk = kargs.get('apk', kargs.get('apk2', 'first'))
        kargs['apk'] = apk
    else:
        apk1 = kargs.get('apk1', 'max')
        kargs['apk'] = apk1

    pk(source, axis=axis, **kargs)

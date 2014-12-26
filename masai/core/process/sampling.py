
#===============================================================================
# undersampling
#===============================================================================
def undersampling(source, N):
    """
    Under sampling the F2 dimension
    """
    #TODO: needs to understand why apk is not perfomring well with undersampling?
    if not source.getIsTransformed():
        source.history.append('undersampling -1')
        nr, nc = source.data.shape
        nc0 = N*int(nc/N) # we need to have nc divide by 2

        source.storePar('TD',nc0)
        source.storePar('SW_h', source.fetchPar('SW_h')/float(N))

        source.data = source.data.iloc[:, :nc0:N]

#===============================================================================
# clip region
#===============================================================================
def clip(source, c0, c1, axis=-1):
    """
    Clip a region of the spectra

    """

    # if not source.getIsTransformed(axis):
    #    return

    source.history.append('clipping %d' % axis)

    if source.is_2D and axis== 0:
        data = source.data.T.copy()
    else:
        data = source.data.copy()
        axis = -1

    # clip the data
    x = data.columns
    select = x[(x>=c0) & (x<=c1)]
    if select.size % 2 != 0: select = select[:-1]
    S = data.ix[:,select]

    if source.getIsTransformed(axis):
        # recompute SI, SW and new apparent SF
        SI0 = data.shape[-1]
        SI = S.shape[-1]
        SW_h = source.fetchPar('SW_h', axis) * float(SI) / float(SI0)
        source.storePar('SW_h', SW_h, axis)
        source.storePar('SI', SI, axis)
        BF1 = source.fetchPar('BF1', axis)
        SFO1 = source.fetchPar('SFO1', axis)
        delta = -SW_h / (float(SI) * BF1)
        SF = SFO1 - (select[-1] + delta * float(SI) / 2.) * BF1 * 1.e-6
        source.storePar('SF', SF, axis)
    else:
        TD = S.shape[-1]
        TDeff = TD
        source.storePar('TD', TD, axis)
        source.storePar('TDeff', TD, axis)

    # save new data
    if source.is_2D and axis == 0:
        source.data = S.T
    else:
        source.data = S

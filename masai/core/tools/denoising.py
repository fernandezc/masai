import numpy as np
from pandas import Series, DataFrame

from core.analysis.pca import PCA


def PCAdenoise(data, quality=99.999999, npc_min=2, verbose=True,
                                                ft=False, axis=None, real=True):

    is_dataframe = False
    if isinstance(data, DataFrame):
        rows = data.index
        cols = data.columns
        data = data.values  # get the ndarray
        is_dataframe = True

    if verbose: print "Denoising using PCA: %s" % ('*' * 59)

    if axis == 0 or axis is None:
        # The data must be transposed for work on axis 0
        S = data.copy().T
        if ft:
            # apply a fft first
            S = np.fft.ifft2(S)
        # Perform PCA
        _pca = PCA(S, treshold=quality, npc_min=npc_min)
        S =_pca.transform()
        if ft:
            # get back the untransformed data (but denoised)
            S = np.fft.fft2(S)
        # untranspose
        S = S.T
        # write some informations
        if verbose:
            print "\tAxis:0",
            if ft:
                print "using fft data,",
            _pca.quality()

    if axis == 1:
        S = data.copy()
    if axis == 1 or axis is None:
        if ft:
            # apply a fft first
            S = np.fft.ifft2(S) 
        _pca = PCA(S, treshold=quality, npc_min=npc_min)
        S = _pca.transform()
        if ft:
            S = np.fft.fft2(S)
        # write some informations
        if verbose:
            print "\tAxis:1",
            if ft:
                print "using fft data,",
            _pca.quality()

    if verbose:
        print '*' * 80

    # return a dataframe when necessary
    if real:
        S = S.real
        print 'real', S[0, 0]
    if is_dataframe:
        S = DataFrame(S, index=rows, columns=cols)
    return S

def svd_denoising(f, orda, rank=0, threshold=1., remove=1, verbose=False):
    """Following the Pascal Man method"""

    si = f.size

    # Create the Hankel matrix (H(i,j) = H(i-1,j+1)
    n1 = si - orda + 1
    if n1 < orda:
        print "orda is too large for this data, orda max : %d" % (si / 2)
        print "corrected to %d"%(si/2-2)
        orda = si/2 -2
        n1 = si - orda + 1

    H = np.zeros((n1, orda), dtype=np.complex128)
    for i in xrange(n1):
        H[i] = f[i:orda + i].copy()

    # SVD decomposition
    U, s, V = np.linalg.svd(H, full_matrices=True)

    # Cleaning and Reconstruction
    S = np.zeros(H.shape, dtype=complex)

    sk = s.copy()
    ev = (sk ** 2)
    evr = ev / ev.sum()
    evrc = evr.cumsum()
    if rank == 0:
        rank = orda
    else:
        sk[rank:] = 0.0
    if threshold < 1.0:
        rank = max(abs(evrc - threshold).argmin() + 1, rank)
        sk[rank:] = 0.0
    if remove == 1:
        power = np.sum((s - sk) ** 2)  # total power removed
        n = len(s) - len(np.nonzero(sk)[0])  # number of entry removed
        if n > 0 :
            sk = np.sqrt((sk ** 2 - power / n).clip(0.0, s[0] ** 2))  # do not forget

    if verbose:
        print "keep rank %d - quality = %.6f %%" % (rank, evrc[ rank - 1] * 100.)

        # sk = np.where(sk < threshold, 0.0, sk)

    S[:orda, :orda] = np.diag(sk)

    nH = np.dot(U, np.dot(S, V))

    # recreate Hankel#     (si / 2 + n + 1, si / 2 - n)
    for i in xrange(orda):
        a = nH[0, i]
        g = 1
        for j in xrange(i):
            a += nH[i - j, j]
            g += 1
        nH[0, i] = a / float(g)

        # last columns
    for i in xrange(1, n1):
        a = nH[i, -1]
        g = 1
        for j in xrange(orda - 1):
            if i + j + 1 > n1 - 1:
                break
            a += nH[i + j + 1, -(2 + j)]
            g += 1
        nH[i, -1] = a / float(g)

    # restaure the denoised signal
    for i in xrange(orda):
        f[i] = nH[0, i]
    for i in xrange(1, n1):
        f[orda + i - 1] = nH[i, -1]

    return f

def SVDdenoise(_data, orda, rank=0, threshold=100., remove=1, n_iter=1,
                                             axis= -1, real=True, verbose=True):
    is_series = False
    is_dataframe = False
    if isinstance(_data, Series):
        index = _data.index
        data = _data.values.copy()  # get the ndarray
        is_series = True

    elif isinstance(_data, DataFrame):
        rows = _data.index
        cols = _data.columns
        data = _data.values.copy()  # get the ndarray
        is_dataframe = True

    else:
        raise Exception("input data for SVDdenoise must be a DataFrame or a series")

    NPT = data.shape[axis]
    if data.ndim == 2:
        nrows, ncols = data.shape

    if verbose:
        print "Denoising using SVD of Hankel matrix: %s" % ('*' * 42)
        print "\tAxis: %d Hankel matrix:[%d,%d] rank=%d or threshold=%.5f" % (
                                  axis, NPT - orda + 1 , orda, rank, threshold)
        print "*"*80

    fullverb = False
    if verbose == 'Full':
        fullverb = True

    threshold = threshold / 100.
    # 1D case
    if len(data.shape) == 1:
        fid = np.fft.ifft(data, NPT)
        for _ in xrange(n_iter):
            fid = svd_denoising(fid, orda, rank, threshold, remove, 
                                                               verbose=fullverb)
        s = np.fft.fft(fid, NPT)
        if real:
            s = s.real
        data = s[:NPT]
        # return a dataframe when necessary
        if is_series:
            data = Series(data, index=index)
        return data

    # 2D case
    if axis == -1:
        axis = 1
    if axis == 0:
        for col in xrange(ncols):
            s = data[:, col]
            fid = np.fft.ifft(s, NPT)
            for _ in xrange(n_iter):
                fid = svd_denoising(fid, orda, rank, threshold, remove, 
                                                               verbose=fullverb)
            s = np.fft.fft(fid)
            if real:
                s = s.real
            data[:, col] = s[:NPT]
    elif axis == 1:
        for row in xrange(nrows):
            s = data[row]
            fid = np.fft.ifft(s, NPT)
            for _ in xrange(n_iter):
                fid = svd_denoising(fid, orda, rank, threshold, remove, 
                                                               verbose=fullverb)
            s = np.fft.fft(fid, NPT)
            if real:
                s = s.real
            data[row] = s[:NPT]

    if is_dataframe:
        data = DataFrame(data, index=rows, columns=cols)

    return data


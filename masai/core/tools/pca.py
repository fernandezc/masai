import numpy as np
from scipy.linalg import svd

from traits.api import HasTraits, Array, Float, Int, Bool

#===============================================================================
# PCA
#===============================================================================
class PCA(HasTraits):

    X = Array
    U = Array
    D = Array
    V = Array
    treshold = Float
    npc_min = Int
    npc_max = Int
    _npc = Int
    center = Bool
    wr = Array  # horizontal axisscale
    t = Array

    def __init__(self, specs, treshold=99.9999, npc_min=2, npc_max=200,
                 center=False):
        """
        initialisation
        """

        self.treshold = treshold / 100.
        self.npc_min = npc_min
        self.npc_max = npc_max
        self.center = center
        if specs is not None:
            self.fit(specs)

    def fit(self, specs):
        """
        Compute PCA
        specs : array
            Contains the array of spectra to decompose
        """
        self.X = specs.copy()

        # centering
        if self.center:
            self.mean_ = np.mean(self.X, axis=0)
            self.X -= self.mean_
        else:
            self.mean_ = 0.

        # SVD decomposition
        self.U, self.D, self.V = svd(self.X, full_matrices=False)
        self.ev = ev = (self.D ** 2) / float(self.X.shape[0])
        self.evr = evr = ev / ev.sum()
        self.evrc = evr.cumsum()

    def reset_npc(self):
        npc = min(abs(self.evrc - self.treshold).argmin() + 1, self.npc_max)
        npc = max(self.npc_min, npc)
        return npc

    def pcamodel(self, npc=None):
        if npc is None:
            npc = self.reset_npc()
        scores = self.U[:, :npc] * self.D[:npc]
        loadings = self.V[:npc, :]
        return (scores, loadings)

    def quality(self):
        npc = self.reset_npc()
        print "with %d PC components - fit = %.6f %%" % (npc, self.evrc[npc - 1] * 100.)
        return npc

    def transform(self, scores=None, loadings=None, npc=None):
        """
        Returned inverse transform data using the given PC's
        """
        if scores is None or loadings is None:
            scores, loadings = self.pcamodel(npc)
        return (np.dot(scores, loadings) + self.mean_)


#===============================================================================
# test
#===============================================================================
if __name__ == '__main__':

    pass


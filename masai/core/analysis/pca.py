import numpy as np
from scipy.linalg import svd

from traits.api import HasTraits, Array, Float, Int, Bool, Str, Property, cached_property

# ===============================================================================
# PCA
# ===============================================================================
class PCA(HasTraits):
    X = Array
    U = Array
    D = Array
    V = Array
    threshold = Float
    npc_min = Int
    npc_max = Int
    npc = Property
    centered = Bool
    message = Str


    def __init__(self, specs,
                 threshold=99.9999, npc_min=2, npc_max=200,
                 centered=False):
        """
        initialisation
        """

        self.threshold = threshold / 100.
        self.npc_min = npc_min
        self.npc_max = npc_max
        self.centered = centered
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
        if self.centered:
            self.mean_ = np.mean(self.X, axis=0)
            self.X -= self.mean_
        else:
            self.mean_ = 0.

        # SVD decomposition
        self.U, self.D, self.V = svd(self.X, full_matrices=False)
        self.ev = ev = (self.D ** 2) / float(self.X.shape[0])
        self.evr = evr = ev / ev.sum()
        self.evrc = evr.cumsum()

    @cached_property
    def _get_npc(self):
        npc = min(abs(self.evrc - self.threshold).argmin() + 1, self.npc_max)
        npc = max(self.npc_min, npc)
        return npc

    def model(self, npc=None):
        if npc is None:
            npc = self.npc
        scores = self.U[:, :npc] * self.D[:npc]
        loadings = self.V[:npc, :]
        return (scores, loadings)

    def quality(self, npc=None):
        if npc is None:
            npc = self.npc
        return self.evrc[npc - 1] * 100.

    def transform(self, scores=None, loadings=None, npc=None):
        """
        Returned inverse transform data using the given PC's
        """
        if npc is None:
            npc = self.npc

        if scores is None or loadings is None:
            scores, loadings = self.model(npc)

        return (np.dot(scores, loadings) + self.mean_)


#===============================================================================
# test
#===============================================================================
if __name__ == '__main__':
    pass


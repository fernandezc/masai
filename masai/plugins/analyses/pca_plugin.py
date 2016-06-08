# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.analyses.pca_plugin
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================
"""
Define pca plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  general imports
# ===============================================================================
# noinspection PyUnresolvedReferences
from traits.api import (HasTraits, Bool, Any, ListStr, List, Either, Dict,
                        Array, Str, Instance, provides, Property)
import traits.has_traits
import pandas as pd
import numpy as np

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.source import Source
from masai.core.utils import htmldoc
from masai.core.analysis.pca import PCA
from masai.core.process import baseline

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Pca'

# ===================================================================================
# Process1d plugin
# ===================================================================================
@provides(IPlugin)
class Pca(HasTraits):
    """ PCA (type: analyses), PCA analysis of a series of spectra.
    """
    # minimal plugin attribute -----------------------------------------------------
    plugin_info = Str("Pca (type: analyses), PCA analysis of a series of spectra.")
    plugin_type = Str("Pca")

    # specific attributes ----------------------------------------------------------
    sources = Either(Source, Dict)

    loadings = Instance(pd.DataFrame)
    scores = Instance(pd.DataFrame)
    message = Str
    npc = Property

    silent = Bool(True)

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    def _repr_html_(self):
        if not self.sources:
            return htmldoc(self.__init__.__doc__)
        else:
            return self.message

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """
        PCA analysis of a list of sources or 2D spectra
        Usage: Pca([[]source[, ]], **kargs)

        Parameters
        ----------
        source: List
            A list of source instances

        Properties
        ----------
        npc: int
            number of components

        Methods
        -------
        quality: usage f.quality([npc])
            return the quality of the decomposition with the given
            number of components
        """
        self.isdic = False
        if args:
            # look in args
            if isinstance(args[0], dict):
                # we assume homogeneous dataframe same index and columns
                self.isdict = True
                self.sources = s = args[0]
                d = []
                self.dindex = dindex = s.keys()
                for k in dindex:
                    ir = s[k].data.values.real
                    m=ir.min()
                    d.append(ir-m)

                self.columns = s[dindex[0]].data.columns
                self.data = data = np.concatenate(d, axis=0)
                parvar = kargs.get('parvar')
                if isinstance(parvar, list) and len(parvar)==len(d):
                    self.index = parvar
                else:
                    self.index = range(len(d))

            elif (isinstance(args[0], Source) and args[0].is_2d):
                self.data = data = args[0].data.values
                self.columns = args[0].data.columns
                self.index = args[0].data.index

        else:
            return

        npc_min = kargs.get('npc_min', 2)
        npc_max = kargs.get('npc_max', 30)
        threshold = kargs.get('threshold', 99.9)
        centered = kargs.get('centered', False)

        self.pca = pca = PCA(data, npc_min=npc_min, npc_max=npc_max,
                             threshold=threshold, centered=centered)

    def model(self, npc=None):
        scores, loadings = self.pca.model(npc=npc)
        # make a dataframe of loadings and scores
        self.loadings = loadings = pd.DataFrame(loadings,
                                                columns=self.columns)
        self.scores = scores = pd.DataFrame(scores) #, index=self.index)
        return loadings, scores

    def transform(self, npc=None):
        if self.loadings is None:
            print('perform pca.model() first!')
            return None
        if npc is None:
            npc = self.npc

        data = self.pca.transform(scores=self.scores.values,
                                  loadings=self.loadings.values, npc=npc)


        if self.isdict:
            newdata = np.split(data, len(self.dindex), axis=0)
            newsources = self.sources.copy()
            for i, dindex in enumerate(self.dindex):
                newsources[dindex].values = newdata[i]
            return newsources
        else:
            return pd.DataFrame(data, index=self.data.index,
                                columns=self.data.columns)

    def quality(self, npc=None):
        if npc is None:
            npc = self.pca.npc
        return self.pca.quality(npc=npc)

    def _get_npc(self):
        return self.pca.npc

    def basecorr(self):
        if self.loadings is None:
            print('perform pca.model() first!')
            return None
        for i in range(self.npc):
            base, w = baseline._polybase(self.loadings.iloc[i].values,
                                         retw=True)
            self.loadings.iloc[i] = self.loadings.values[i] - base

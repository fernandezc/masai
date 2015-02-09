# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
# ===============================================================================

from __future__ import print_function
import unittest
import os

from masai.api import Bruker, Process, Fit


def read_fid(DATADIR, name, expno):
    path = os.path.join(DATADIR, 'user', 'nmr', name, str(expno))
    source = Bruker(path)
    return source


def read_series(DATADIR, name, expnos, varpars):
    return Bruker(data_dir=DATADIR,
                  user='user', name=name, expnos=expnos, varpars=varpars)


class FittingTests(unittest.TestCase):
    def setUp(self):
        # we need to find where this file is located
        here = (__file__.split(os.path.sep))[:-2]
        here = "/" + os.path.join(*here)
        self.DATADIR = os.path.expanduser(os.path.join(here, 'doc', 'tutorial', 'exemples'))

    def test_fitting(self):
        source = read_fid(self.DATADIR, 'SPEC1D', 10)
        Process(source, transform="em 50; zf 1; ft --tdeff 512")
        Process(source, transform="pk --auto --fit_phc1 --bound_phc1=100; ab --mode linear")

        script = """
        COMMON:
        # common parameters ex.
        # $ gwidth: 1.0, 0.0, none
          $ gratio: 0.5, 0.0, 1.0

        MODEL: LINE_1
        shape: voigtmodel
            * ampl:  1.0, 0.0, none
            $ pos:   4.5, -100.0, 100.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_2
        shape: voigtmodel
            $ ampl:  0.2, 0.0, none
            $ pos:   1.6, 0, 2.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_3
        shape: voigtmodel
            $ ampl:  0.18, 0.0, none
            $ pos:   7.9, 5.0, 10.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_4
        shape: voigtmodel
            $ ampl:  0.1, 0.0, none
            $ pos:   14, 10.0, 20.0
            > ratio: gratio
            $ width: 6.0, 0, 20

        MODEL: LINE_5
        shape: voigtmodel
            $ ampl:  0.1, 0.0, none
            $ pos:   -5, -8.0, -1.0
            > ratio: gratio
            $ width: 6.0, 0, 10
        """
        f = Fit(source, script)
        f.run(maxiter=1000, every=200)
        
if __name__ == '__main__':
    unittest.main()

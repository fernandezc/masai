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

import masai

from masai.api import Bruker, Process, Multiplex


def read_fid(DATADIR, name, expno):
    path = os.path.join(DATADIR, 'user', 'nmr', name, str(expno))
    source = Bruker(path)
    return source


def read_series(DATADIR, name, expnos, varpars):
    return Bruker(data_dir=DATADIR,
                  user='user', name=name, expnos=expnos, varpars=varpars)


class BasicTests(unittest.TestCase):
    def setUp(self):
        # we need to find where this file is located
        here = (__file__.split(os.path.sep))[:-2]
        here = "/" + os.path.join(*here)
        self.DATADIR = os.path.expanduser(os.path.join(here, 'doc', 'tutorial', 'exemples'))

    def test_bruker_plugin(self):
        f1d = read_fid(self.DATADIR, 'SPEC1D', 10)
        f2d = read_fid(self.DATADIR, 'HMQC', 31)
        f3d = read_series(self.DATADIR, 'CP', range(80, 95), "P 15")
        row = f2d.get_row(0, byindex=False)
        col = f2d.get_col(0, width=.1, byindex=False)
        row = f2d.get_row(1, byindex=True)
        col = f2d.get_col(1, byindex=True)

        user = 'user'
        name = 'MULTIPLEX'
        expno = 4
        path = os.path.join(self.DATADIR, user, 'nmr', name, str(expno))
        source = Multiplex(path)
        # possibly correct the F1 encoding
        source.par2.encoding = 'STATES'

    def test_process_plugin(self):
        f1d = read_fid(self.DATADIR, 'SPEC1D', 10)
        Process(f1d, transform="em 100; zf 0; ft")
        Process(f1d, transform="pk --fit_phc1")

        f2d = read_fid(self.DATADIR, 'HMQC', 31)
        d6 = f2d.par.D[6] * 1.e3  # in ms
        Process(f2d, transform='em 50 --shifted %.3f' % d6)
        Process(f2d, transform='zf 0')
        Process(f2d, transform='ft')
        Process(f2d, transform='pk 0 0 --auto --fit_phc1 --pivot 38.4 --shifted %.3f' % d6)
        Process(f2d, axis=0, transform='em 20')
        Process(f2d, axis=0, transform='zf 0')
        Process(f2d, axis=0, transform='ft')
        Process(f2d, axis=0, transform='pk 180 0 --auto --fit_phc1 --select cols')

    def test_STMAS(self):
        source = read_fid(self.DATADIR, 'STMAS', '15')

        Process(source, transform='em 10 --k_shifted 1 2 --verbose')
        Process(source, transform='us 2')
        Process(source, transform='zf 1')
        Process(source, transform='ft')
        Process(source, transform='cl -200 0')
        Process(source, transform='pk -256.364 0 --fit_phc1  --auto')
        Process(source, "ab --mode poly")

        Process(source, axis=0, transform='lp 128 --tdeff 108')  # must be before zf in any case
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 3')
        Process(source, axis=0, transform='ft --shear 1 2 --shift 0')
        Process(source, axis=0, transform='pk 0 0 --auto --fit_phc1 --select cols --verbose')

        #TODO make tests also for plotting


if __name__ == '__main__':
    pass

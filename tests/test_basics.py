# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
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

from masai.api import *

def read_fid(DATADIR, name, expno):
    path = os.path.join(DATADIR, 'user', 'nmr', name, str(expno))
    source = Bruker(path)
    return source


class BasicTests(unittest.TestCase):
    def setUp(self):
        self.DATADIR = os.path.expanduser(os.path.join('..', 'doc', 'tutorial', 'exemples'))

    def test_bruker_plugin(self):
        f1d = read_fid(self.DATADIR, 'SPEC1D', 10)
        f2d = read_fid(self.DATADIR, 'HMQC', 31)

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

        # Process(source, axis=0, transform='lp 128 --tdeff 108') # must be before zf in any case
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 3')
        Process(source, axis=0, transform='ft --shear 1 2 --shift 0')
        Process(source, axis=0, transform='pk 0 0 --auto --fit_phc1 --select cols --verbose')

        # """
        # test plugin
        # """
        # DATADIR = '../../../doc/tutorial/exemples'
        # sources = Bruker(data_dir=DATADIR,
        #         user = 'user',
        #         name = 'CP',
        #         expnos = range(80,95,1),
        #         varpars = "P 15")
        #
        # #print(sources)
        # import sys
        # sys.exit()
        #
        # # need this to test  directly
        # import pylab as pl
        # from masai.plugins.catalog import Process, Plot2d
        #
        # Bruker = Bruker()
        # Process = Process()
        # Plot2d = Plot2d()
        #
        # path = '/Users/christian/Dropbox/PycharmProjects/masai/exemples/user/nmr/HMQC/31'
        # fd = Bruker(path)
        # print(fd)
        # d6 = fd.par.D[6]*1e3
        # Process(fd, transform='em 5000 --shifted %.2f'%d6)
        # #Process(fd, transform='zf 0; ft; pk -78 -10')
        # Plot2d(fd)
        #
        #
        # import sys
        # sys.exit()
        #
        # DATADIR = os.path.expanduser(os.path.join('~','Dropbox','storage','bruker','data'))
        #
        # user = 'christian'
        # name = 'example1'
        # expno = 10
        # path = os.path.join(DATADIR,user,'nmr',name, str(expno))
        # fd = Bruker(path, dgfilter=True)
        #
        # user = 'aurelie'
        # name = 'RbNO3'
        # expno = 15
        # path = os.path.join(DATADIR,user,'nmr',name, str(expno))
        # fd = Bruker(path, dgfilter=True)
        #
        #
        #
        # #_----------------------------------
        # import warnings
        # warnings.filterwarnings("ignore")
        #
        # Bruker = Bruker()
        #
        #
        #
        # sources = Bruker(
        #             data_dir=DATADIR,
        #             user = 'siblanixenon',
        #             name = 'NaY1-QIN',
        #             expnos = [25,49,65,85,105,124,144,164,184,194,209,224,242],
        #             varpars = "TE",
        #             )
        #
        # sources = Bruker(
        #              data_dir=DATADIR,
        #              user = 'siblanixenon',
        #              name = 'NaY1-QIN',
        #              expnos = [25,49,65,85,105,124,144,164,184,194,209,224,242],
        #              varpars = "TE",
        #              )
        #
        # sources.data.T.plot(legend=False)
        # #print(sources)
        # pl.show()
        #
        #
        # import sys
        # sys.exit()
        #
        # user = 'christian'
        # name = 'example1'
        # expno = 10
        # path = os.path.join(DATADIR,user,'nmr',name, str(expno))
        # fd = Bruker(path, dgfilter=True)
        # fd.data.T.plot(legend=False)
        #
        # user = 'aurelie'
        # name = 'RbNO3'
        # expno = 15
        # path = os.path.join(DATADIR,user,'nmr',name, str(expno))
        # fd = Bruker(path, dgfilter=True)
        #
        # row = fd.row(0, byindex=False)
        # row.T.plot(legend=False)
        #
        # col = fd.col(0, width=.1, byindex=False)
        # col.T.plot(legend=False)
        #
        # pl.show()
        #
        # #TODO make tests also for plotting


if __name__ == '__main__':
    unittest.main()

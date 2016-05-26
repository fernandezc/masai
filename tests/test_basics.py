# -*- coding: utf-8 -*-

"""

Test basics

"""

# Copyright (c) C. Fernandez' group @ LCS.
# Distributed under the terms of the CeCILL-B License.

from __future__ import print_function

# -----------------------------------------------------------------------------
# Globals
# -----------------------------------------------------------------------------

# stdlib
import unittest
import os

# third-party

# -----------------------------------------------------------------------------
# Local utilities and constants
# -----------------------------------------------------------------------------
from masai.api import *


def read_fid(DATADIR, name, expno, sr2=None, sr1=None):
    path = os.path.join(DATADIR, 'user', 'nmr', name, str(expno))
    source = Bruker(path, sr2 = sr2, sr1= sr1)
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

        SR = 1551.
        source = read_fid(self.DATADIR, 'EXSY', 10, sr2 = SR, sr1= SR)
        Process(source, transform='em 1 --k_shifted 1 0')
        Process(source, transform='zf 0')
        Process(source, transform='ft')
        Process(source, transform='cl -4 13')
        ph = 80
        Process(source, transform='pk %.3f 0 --pivot 6.85' % ph)
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 2')
        Process(source, axis=0, transform='ft')
        Process(source, axis=0, transform='pk 2.5 0 --select cols --pivot 0')
        Process(source, axis=0, transform='cl -4 13')

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


    def test_tutorial_basics(self):

        # get the doc for a plugin
        Bruker()

        # specifications for the data to open
        DATADIR = os.path.join(MASAI,'doc', 'tutorial', 'exemples')
        user = 'user'
        name = 'SPEC1D'
        expno = 10

        # load it
        path = os.path.join(DATADIR, user,'nmr',name, str(expno))
        source = Bruker(path)
        source = Bruker(data_dir=DATADIR, user = 'user', name = 'SPEC1D', expno = 10)
        source = Bruker(data_dir=DATADIR, user = 'user', name = 'SPEC1D', expno = 10)
        print('All P:', source.par.P)
        print('P2:', source.par.P[2])
        print('DE:', source.par.DE)
        p = Plot1d(source, display=False)  # basic default plot

        p = Plot1d(source, imag=True,
               figsize=(3.0,2.0), dpi=100, fontsize=8, xlim=(0, 4), c='brown',
               display=False)

        p = Plot1d(source, hold=True, color='brown',
                   display=False)
        p = Plot1d(source, axe=p.axe, imag=True, xlim=(0,4), color='orange', ls='--',
                   display=False)
                        # we pass a reference to the previous axe (a)

        p = Plot1d(source, color='orange', ls='-', lw='5', alpha='.4', xlim=(0,4),
                   commands = ["text(2.0, -0.1, 'xxxxx', ha='center', va='center', color='b' )",
                               "text(2.0, -0.1, 'xxxxx', ha='center', va='center', color='g', rotation=90 )",
                              ],
                   display=False)

        p = Plot1d(source, xlabel='decay (ms)', ylabel='a.u', xlim=(0,10),
                   display=False)

        source = Bruker(path)
        Process(source, transform='em 50; zf 0; ft; pk -178 -10')
        p = Plot1d(source, display=False)

        source = Bruker(path)  # reset to original data
        Process(source, transform="em 50; zf 1; ft") # 1 zerofilling
        Process(source, "cl -50 50")
        Process(source, "pk --auto")
        ax = Plot1d(source, display=False)

        source = Bruker(path)  # reset to original data
        Process(source, transform="em 50; zf 1; ft") # 1 zerofilling
        Process(source, "cl -50 50")
        Process(source, "pk --auto --fit_phc1 --bound_phc1 180 --pivot 4.5")
        a = Plot1d(source, display=False) #, xlim=(30,-10))

        p = Plot1d(source, hold=True, display=False)
        Process(source, "ab --mode poly")
        Plot1d(source, axe=p.axe, display=False)

        source = Bruker(path)  # reset to original data
        Process(source, transform="em 50; zf 1; ft; pk --auto --fit_phc1 --bound_phc1 360 --pivot 4.5")
        p = Plot1d(source, hold=True, display=False)
        Process(source, "ab --mode poly")
        p = Plot1d(source.par.baseline, axe=p.axe, xlim=(30,-10), reverse=True, display=False)

if __name__ == '__main__':
    pass

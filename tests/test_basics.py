from __future__ import print_function
import unittest
import os

from masai.api import *

def read_fid(DATADIR, name,expno):
    path = os.path.join(DATADIR,'user','nmr',name, str(expno))
    source = Bruker(path)
    return source

class BasicTests(unittest.TestCase):

    def setUp(self):
        self.DATADIR = os.path.expanduser(os.path.join('..','doc','tutorial','exemples'))

    def test_1D(self):
        f1d = read_fid(self.DATADIR, 'SPEC1D', 10)
        Process(f1d, transform="em 100; zf 0; ft")
        Process(f1d, transform="pk --fit_phc1")

    def test_2D(self):
        f2d = read_fid(self.DATADIR, 'HMQC', 31)
        d6 = f2d.par.D[6]*1.e3 # in ms
        Process(f2d, transform = 'em 50 --shifted %.3f'%d6)
        Process(f2d, transform = 'zf 0')
        Process(f2d, transform = 'ft')
        Process(f2d, transform = 'pk 0 0 --auto --fit_phc1 --pivot 38.4 --shifted %.3f'%d6)
        Process(f2d, axis=0, transform='em 20')
        Process(f2d, axis=0, transform='zf 0')
        Process(f2d, axis=0, transform='ft')
        Process(f2d, axis=0, transform='pk 180 0 --auto --fit_phc1 --select cols')

    def test_STMAS(self):
        source = read_fid(self.DATADIR, 'STMAS', '15')

        Process(source, transform = 'em 10 --k_shifted 1 2 --verbose')
        Process(source, transform = 'us 2')
        Process(source, transform = 'zf 1')
        Process(source, transform = 'ft')
        Process(source, transform = 'cl -200 0')
        Process(source, transform = 'pk -256.364 0 --fit_phc1  --auto')

        #Process(source, axis=0, transform='lp 128 --tdeff 108') # must be before zf in any case
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 3')
        Process(source, axis=0, transform='ft --shear 1 2 --shift 0')
        Process(source, axis=0, transform='pk 0 0 --auto --fit_phc1 --select cols --verbose')

    #TODO make tests also for plotting

if __name__ == '__main__':
    unittest.main()

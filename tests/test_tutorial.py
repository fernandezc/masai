# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.tests.test_tutorial.py
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
"""
Here we use runipy to generate all output from the tutorial notebooks
and may be detect exceptions
"""
from __future__ import print_function

import unittest
import os
import codecs

from IPython import nbformat
from IPython.nbconvert.preprocessors.execute import ExecutePreprocessor

class TutorialTests(unittest.TestCase):

    def setUp(self):
        masai_dir = os.environ.get('MASAI')
        assert os.path.exists(masai_dir)
        self.path = os.path.join(masai_dir, 'doc', 'tutorial')
        self.pp = ExecutePreprocessor()
        self.pp.timeout = 300  # seconds
        self.pp.interrupt_on_timeout = True

    def test_ipynb(self):

        ldir = os.listdir(self.path)
        for filename in ldir:
            if os.path.splitext(filename)[1]==".ipynb" \
                    and not filename.startswith('~'):
                print("\n >>>>>> running %s" % filename)
            else:
                continue
            ipynb = os.path.join(self.path, filename)
            nb = nbformat.read(open(ipynb), as_version=4)
            nb_executed, resources = self.pp.preprocess(nb, resources={})
            for cell in nb_executed.cells:
                # print(cell)
                if "outputs" in cell.keys():
                    outs = cell["outputs"]
                    for out in outs:
                        if "output_type" in out and \
                                    str(out["output_type"])=="error":
                            print(cell)
                            print(outs)
                            print(out['traceback'])
                            raise Exception(out['ename'])


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

from runipy.notebook_runner import NotebookRunner
from IPython.nbformat.current import read, write

class TutorialTests(unittest.TestCase):

    def setUp(self):
        masai_dir = os.environ.get('MASAI')
        assert os.path.exists(masai_dir)
        self.path = os.path.join(masai_dir, 'doc', 'tutorial')

    def test_ipynb(self):
        ldir = os.listdir(self.path)
        nb=0
        for item in ldir:
            if os.path.splitext(item)[1]==".ipynb":
                ipynb = os.path.join(self.path, item)
                print("\n >>>>>> running %s" % ipynb)
                nb += 1
            else:
                continue
            notebook = read(open(ipynb), 'json')
            r = NotebookRunner(notebook)
            r.run_notebook()
            write(r.nb, open(ipynb, 'w'), 'json')

        # is there any notebook tested?
        assert nb>0


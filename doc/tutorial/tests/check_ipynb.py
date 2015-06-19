# -*- coding: utf-8 -*-

# =======================
#  Test ipython notebook
# =======================

# based on recipes found on:
# https://github.com/maxalbert/auto-exec-notebook

from __future__ import print_function

import os
import codecs

from IPython import nbformat
from IPython.nbconvert.preprocessors.execute import ExecutePreprocessor

pp = ExecutePreprocessor()
pp.timeout = 300  # seconds
pp.interrupt_on_timeout = True

def check_ipynb(path):
    ldir = os.listdir(path)
    for filename in ldir:
        if os.path.splitext(filename)[1]==".ipynb" \
                and not filename.startswith('~'):
            print("\n >>>>>> running %s" % filename)
        else:
            continue
        ipynb = os.path.join(path, filename)
        nb = nbformat.read(open(ipynb), as_version=4)
        nb_executed, resources = pp.preprocess(nb, resources={})
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


masai_dir = os.environ.get('MASAI')
assert os.path.exists(masai_dir)

path = os.path.join(masai_dir, 'doc', 'tutorial')
check_ipynb(path)

#!/usr/bin/env python
"""strip outputs from an IPython Notebook
Opens a notebook, strips its output, and writes the outputless version to the original file.
Useful mainly as a git pre-commit hook for users who don't want to track output in VCS.
This does mostly the same thing as the `Clear All Output` command in the notebook UI.
Adapted from rom https://gist.github.com/minrk/6176788 to work with
git filter driver

Add this to your .git/config:

[filter "stripoutput"]
    clean = "/path/to/nbconvert/nbstripout"
and a attributes file with

*.ipynb filter=stripoutput


to put in .git/info

#You may need to do this for your script to work with GitX or Tower:
#sys.path.append("/Users/chris/anaconda/envs/conda/lib/python2.7/site-packages")

"""
import sys

from IPython.nbformat import current

def strip_output(nb):
    """strip the outputs from a notebook object"""
    for cell in nb.worksheets[0].cells:
        if 'outputs' in cell:
            cell['outputs'] = []
        if 'prompt_number' in cell:
            cell['prompt_number'] = ""
    return nb

if __name__ == '__main__':
    nb = current.read(sys.stdin, 'json')
    nb = strip_output(nb)
    current.write(nb, sys.stdout, 'json')

    # import io
    # import os
    #
    # filename = sys.argv[1]
    # with io.open(filename, 'r', encoding='utf8') as f:
    #     nb = current.read(f, 'json')
    # nb = strip_output(nb)
    # with io.open(filename, 'w', encoding='utf8') as f:
    #     current.write(nb, f, 'json')
    # path, ext = os.path.splitext(filename)
    # filename_out = os.path.join(path,'_clean',ext)

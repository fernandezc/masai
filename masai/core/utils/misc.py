# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.utils.misc
# ===============================================================================
# Copyright (C) 2012-2013 Christian Fernandez
#       Laboratoire Catalyse et Spectrochimie, Caen, France.
#       christian.fernandez@ensicaen.fr
# This software is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# You can  use, modify and/ or redistribute the software under
# the terms of the CeCILL-B license as circulated by CEA, CNRS and INRIA
# at the following URL "http://www.cecill.info".
# See Licence.txt in the main masai source directory
#===============================================================================
""".. _misc:

Various methods used in other part of the program

"""
import re

import pandas as pd

# import string

def position2index(data, pos):
    """
    Return the index of a position along the last axis (columns) of a DataFrame

    Parameters
    ----------
        data: DataFrame
        pos: float
            a column position in the column units

    Returns
    -------
        index: int
            the index of the position

    """
    if isinstance(data, pd.DataFrame):
        a = data.columns
    else:
        a = data
    if a[0] < a[-1]:
        id = a[a <= pos].size
    else:
        id = a[a >= pos].size

    # check closest
    dpos0 = 1.e15
    idx = 0
    id0 = id
    for i in [-1, 0, +1]:
        id = max(min(id0 + i, a.size - 1), 0)
        dpos = abs(pos - a[id])
        if dpos < dpos0:
            idx = id
            dpos0 = dpos
    return idx


##### Formating #####

def htmldoc(text):
    p = re.compile("^(?P<name>.*:)(.*)", re.MULTILINE)  #To get the keywords
    html = p.sub(r'<b>\1</b>\2', text)
    html = html.replace('-', '')
    html = html.split('\n')
    while html[0].strip() == '':
        html = html[1:]  # suppress initial blank lines

    for i in range(len(html)):
        html[i] = html[i].strip()
        if i == 0:
            html[i] = "<h3>%s</h3>" % html[i]
        html[i] = html[i].replace('Parameters', '<h4>Parameters</h4>')
        html[i] = html[i].replace('Properties', '<h4>Properties</h4>')
        html[i] = html[i].replace('Methods', '<h4>Methods</h4>')
        if html[i] != '':
            if "</h" not in html[i]:
                html[i] = html[i] + '<br/>'
            if not html[i].strip().startswith('<'):
                html[i] = "&nbsp;&nbsp;&nbsp;&nbsp;" + html[i]
    html = "".join(html)

    return html


if __name__ == '__main__':
    text = """
    Fit a 1D, 2D or a list of sources

        Usage: Fit([[]source[, ]], **kargs)

        Parameters
        ----------
        source: Source instance
            The data to fit (may be a list of sources)
        mode: Str, optional
            Reserved - not used for now
        kind: Str, optional
            Contain information about the specific kind of experiment
            for now, only 'DQ-SQ' is implemented

        Properties
        ----------
        fp: Dict
            Fit parameters dictionary (read-only, but individual elements of the dict can be changed)
        script: Str
            A string representation of the fp dict, which can be used as input for other a fit

        Methods
        -------
        run: Method
            Run a fit, use run() without arguments to get help on this method

        """

    html = htmldoc(text)
    print(html)


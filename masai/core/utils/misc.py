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


def is_sequence(arg):
    """
    Determine if an object is iterable but is not a string
    """
    return not hasattr(arg, 'strip') and  hasattr(arg, "__iter__")

def srepr(arg):
    if is_sequence(arg):
        return '<' + ", ".join(srepr(x) for x in arg) + '>'
    return repr(arg)


def makestr(l):
    """
    make a string from a list of string
    """

    if is_sequence(l):
        l = " ".join(map(str,l))
    l = l.replace('$','')
    l = l.replace(' ','\ ')
    l = r'$%s$' % l
    return l

if __name__ == '__main__':


    # Test html_doc
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
    assert(html=="<h3>Fit a 1D, 2D or a list of sources</h3><b>        Usage:</b> Fit([[]source[, ]], **kargs)<br/><h4>Parameters</h4><b>        source:</b> Source instance<br/>&nbsp;&nbsp;&nbsp;&nbsp;The data to fit (may be a list of sources)<br/><b>        mode:</b> Str, optional<br/>&nbsp;&nbsp;&nbsp;&nbsp;Reserved  not used for now<br/><b>        kind:</b> Str, optional<br/>&nbsp;&nbsp;&nbsp;&nbsp;Contain information about the specific kind of experiment<br/>&nbsp;&nbsp;&nbsp;&nbsp;for now, only 'DQSQ' is implemented<br/><h4>Properties</h4><b>        fp:</b> Dict<br/>&nbsp;&nbsp;&nbsp;&nbsp;Fit parameters dictionary (readonly, but individual elements of the dict can be changed)<br/><b>        script:</b> Str<br/>&nbsp;&nbsp;&nbsp;&nbsp;A string representation of the fp dict, which can be used as input for other a fit<br/><h4>Methods</h4><b>        run:</b> Method<br/>&nbsp;&nbsp;&nbsp;&nbsp;Run a fit, use run() without arguments to get help on this method<br/>")


# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.processes.process_plugin
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
"""
Define process plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  general imports
# ===============================================================================
# noinspection PyUnresolvedReferences
from traits.api import (HasTraits, Bool, Any, Str, provides)
import traits.has_traits
import re

# ===============================================================================
# local imports
# ===============================================================================
from masai.plugins import IPlugin
from masai.core.api import *

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Process'

all = ['Process', ]

# ===================================================================================
# Process plugin
# ===================================================================================
@provides(IPlugin)
class Process(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Process (type: processes), output a transformed 1D source.")
    type = Str("processes")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def print_doc(self):
        return self.description

    # -------------------------------------------------------------------------------
    def help(self, command=None):
        from IPython.display import HTML
        if not command:
            return HTML(self._repr_html_())
        else:
            try:
                command(options="--help")
            except TypeError:
                try:
                    exec("%s(options='--help')" % command)
                except SystemExit:
                    pass

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """Process a 1D source

        Parameters
        ----------
        fd: source (such as NMR or IR source)
        transforms: str
            line(s) of commands among: em, zf, ft, pk, ab
            type Process.<command_name> to get help on individual command
        """
        #todo: add new commands
        # object to apply the transform
        if args:
            # look in args
            source = args[0]
        else:
            raise ValueError("no valid source specified\nUsage: Process(fd)")

        # type of transform, by default em+zf ft+apk
        if len(args)>1:
            transforms = args[1]
        else:
            transforms = kargs.get('transform', '-h; em 10; zf; ft; pk --auto')

        rlines = re.compile(r"[ \t\n\r]+")
        transforms = rlines.sub(" ", transforms)
        transforms = transforms.lower().split(';')

        axis = kargs.get('axis', -1 )

        for trans in transforms:

            p = re.compile('^(\w*)\s(.*)')
            m = p.match(trans.strip())
            if m :
                command = m.group(1)
                options = m.group(2)
            else:
                command = trans
                options = ''

            if 'em' in command:
                em(source, options=options, axis=axis)

            elif 'zf' in command:
                zf(source, options=options, axis=axis)

            elif 'cs' in command:
                cs(source, options=options, axis=axis)

            elif 'ls' in command:
                ls(source, options=options, axis=axis)

            elif 'rs' in command:
                rs(source, options=options, axis=axis)

            elif 'us' in command:
                us(source, options=options, axis=axis)

            elif 'cl' in command:
                cl(source, options=options, axis=axis)

            elif 'ft' in command:
                ft(source, options=options, axis=axis)

            elif 'pk' in command:
                pk(source, options=options, axis=axis)

            elif 'ab' in command:
                ab(source, options=options, axis=axis)

    def _repr_html_(self):
        import re
        import string
        p = re.compile("(.*:)?(.*)", re.MULTILINE)
        m = p.findall(self.__call__.__doc__)
        html = ""
        start = False
        for mg in m:
            if mg[0].strip() != "":
                html += "<b>%s</b>:%s<br/>" % mg
            else:
                s = string.join(mg," ")
                if s.strip() and start:
                    html += "&nbsp;"*4+"%s<br/>" % s
                if "----" in s:
                    start = True
                    html = "<h3>Parameters</h3>"
        return html


if __name__ == '__main__':
    import sys, os
    from masai.plugins.list import Bruker
    from masai.plugins.list import Plot1d, Plot2d

    Bruker = Bruker()
    Plot1d = Plot1d()
    Plot2d = Plot2d()
    Process = Process()

    def test_2Dprocess(path2D):

        source = Bruker(path2D)

        d6 = source.par.D[6]*1.e3 # in ms

        Process(source, transform = 'em 50 --shifted %.3f'%d6)
        Process(source, transform = 'zf 0')
        Process(source, transform = 'ft')
        Process(source, transform = 'pk 0 0 --auto --fit_phc1 --pivot 38.4 --shifted %.3f'%d6)

        #p_2 = Plot2d(source, figsize=(3,2.7), ylim=(0, 5))

        row = source.row(0, byindex=True)
        #p_1 = Plot1d(row)

        Process(source, axis=0, transform='em 0')
        Process(source, axis=0, transform='zf 0')
        Process(source, axis=0, transform='ft')
        Process(source, axis=0, transform='pk 180 0 --auto --select cols --verbose')

        #p_2 = Plot2d(source, figsize=(3,2.7))

        row = source.row(20.87,width=1.)
        p_1 = Plot1d(row)

    #tests
    here = os.path.split(__file__)
    DATADIR = os.path.join(here[0],'..','..','..','exemples')
    print (DATADIR)
    user = 'user'
    name = 'HMQC'
    expno = 31
    path2D = os.path.join(DATADIR, user,'nmr',name, str(expno))

    test_2Dprocess(path2D)


    sys.exit()

    # global
    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd = Bruker(path)
    Process(fd, 'us 1')
    Process(fd, "em 50")
   # Process(fd, transform = 'zf 0; cs %.3f --time'%0)
   # Process(fd, "ft")
   # Process(fd, 'cl -60 60')
    #Process(fd, "pk -25 0 --auto --fit_phc1 --ediff 2 --mode negmin+entropy --gamma 0.0001 --bound_phc1 200 --verbose")
   # Process(fd, "ab --mode poly")
    Plot1d(fd)
    Process(fd, transform = 'cs -%.3f --time'%10.01)
    Plot1d(fd)
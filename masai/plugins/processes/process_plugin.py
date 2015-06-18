# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.processes.process_plugin
# ===============================================================================
# Copyright (C) 2015 Christian Fernandez
# Laboratoire Catalyse et Spectrochimie, Caen, France.
# christian.fernandez@ensicaen.fr
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
from masai.core.utils import htmldoc

traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Process'

# ===================================================================================
# Process plugin
# ===================================================================================
@provides(IPlugin)
class Process(HasTraits):
    # minimal plugin attribute -----------------------------------------------------
    plugin_info = Str("Process (type: processes), output a transformed 1D source.")
    plugin_type = Str("processes")

    # specific attributes ----------------------------------------------------------
    nooutput = Bool
    silent = Bool(True)

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # *******************************************************************************
    # Private functions
    # *******************************************************************************
    def _repr_html_(self):
        if self.nooutput:
            return ""

        if not self.source:
            return htmldoc(self.__init__.__doc__)
        else:
            return self.message

    def __repr__(self):
        return self.__init__.__doc__

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """
        Process
        Process a source

        Parameters
        ----------
        fd: source (such as NMR or IR source)
        transforms: str
            line(s) of commands among: em, zf, ft, pk, ab
            type Process(<command_name>) to get help on individual command
        """
        #todo: add new commands
        #todo: solve help problems

        # object to apply the transform

        if args:
            # look in args
            self.source = source = args[0]

            if not isinstance(source, IPlugin):
                command = source
                try:
                    command(options="--help")
                except TypeError:
                    try:
                        exec ("%s(None, options='-h')" % command)
                    except SystemExit:
                        pass
                self.nooutput = True
                return

        else:
            return

        self.nooutput = True

        # type of transform, by default em+zf ft+apk
        if len(args) > 1:
            transforms = args[1]
        else:
            transforms = kargs.get('transform', '-h; em 10; zf; ft; pk --auto')

        rlines = re.compile(r"[ \t\n\r]+")
        transforms = rlines.sub(" ", transforms)
        transforms = transforms.lower().split(';')

        axis = kargs.get('axis', -1)

        for trans in transforms:

            p = re.compile('^(\w*)\s(.*)')
            m = p.match(trans.strip())
            if m:
                command = m.group(1)
                options = m.group(2)
            else:
                command = trans
                options = ''

            if 'em' in command:
                em(source, options=options, axis=axis)

            if 'bc' in command:
                bc(source, options=options, axis=axis)

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

            elif 'pc' in command:
                pc(source, options=options, axis=axis)

            elif 'mc' in command:
                mc(source, options=options, axis=axis)

            elif 'lp' in command:
                lp(source, options=options, axis=axis)


if __name__ == '__main__':
    import sys, os
    from masai.api import Bruker
    from masai.api import Plot1d, Plot2d

    Process("ft")

    sys.exit()


    def test_2D_process(path):
        source = Bruker(path)

        d6 = source.par.D[6] * 1.e3  # in ms

        Process(source, transform='em 50 --shifted %.3f' % d6)
        Process(source, transform='zf 0')
        Process(source, transform='ft')
        Process(source, transform='pk 0 0 --auto --fit_phc1 --pivot 38.4 --shifted %.3f' % d6)

        #p_2 = Plot2d(source, figsize=(3,2.7), ylim=(0, 5))

        row = source.row(0, byindex=True)
        #p_1 = Plot1d(row)

        Process(source, axis=0, transform='em 0')
        Process(source, axis=0, transform='zf 0')
        Process(source, axis=0, transform='ft')
        Process(source, axis=0, transform='pk 180 0 --auto --select cols --verbose')

        p_2 = Plot2d(source, figsize=(3, 2.7))

        row = source.row(-27, width=1.)
        p_1 = Plot1d(row)

    def test_STMAS_process(path):
        source = Bruker(path)

        Process(source, transform='em 10 --k_shifted 1 2 --verbose')
        Process(source, transform='us 2')
        Process(source, transform='zf 1')
        Process(source, transform='ft')
        Process(source, transform='cl -200 0')
        Process(source, transform='pk -256.364 0 --fit_phc1  --auto')

        #Process(source, axis=0, transform='lp 128 --tdeff 108') # must be before zf in any case
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 3')
        Process(source, axis=0, transform='ft --shear 1 2 --shift 0')
        Process(source, axis=0, transform='pk 0 0 --auto --fit_phc1 --select cols --verbose')

        p_2 = Plot2d(source, proj=False, figsize=(3, 2.7))  #, xlim=(-20,-50))

    #tests
    here = os.path.split(__file__)
    DATADIR = os.path.join(here[0], '..', '..', '..', 'exemples')
    print(DATADIR)
    user = 'user'
    name = 'HMQC'
    expno = 31
    path2D = os.path.join(DATADIR, user, 'nmr', name, str(expno))

    test_2D_process(path2D)

    name = 'STMAS'
    expno = 15
    pathSTMAS = os.path.join(DATADIR, user, 'nmr', name, str(expno))

    test_STMAS_process(pathSTMAS)

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
    Process(fd, transform='cs -%.3f --time' % 10.01)
    Plot1d(fd)

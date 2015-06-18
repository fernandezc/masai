#! /usr/bin/env python
# ===============================================================================
# modelling.sequence
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

""" docstring """
import os

from traits.api import (HasTraits, Constant,
                        List, Code, Button, Str)

try:  #hack to avoid error with py2app applications
    _scriptdir = os.path.join(os.path.dirname(__file__), 'script')
    _pulse_sequence_list = os.listdir(_scriptdir)
except OSError:
    _scriptdir = 'script'
    _pulse_sequence_list = os.listdir(_scriptdir)


class SequenceAvailable(HasTraits):
    type = List(_pulse_sequence_list)


sequence_available = SequenceAvailable()


class PulseSequence(HasTraits):
    type = Str('ideal_pulse')
    available = Constant(sequence_available)
    _name = Str()

    simpson_script = Code()
    script = Code()

    def _script_default(self):
        with open(os.path.join(_scriptdir, self.type), 'rb') as f:
            return f.read()

    def _type_changed(self):
        with open(os.path.join(_scriptdir, self.type), 'rb') as f:
            self.script = f.read()


#===============================================================================
if __name__ == '__main__':
    PulseSequence().configure_traits()

    pass

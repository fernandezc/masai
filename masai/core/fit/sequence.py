#! /usr/bin/env python
#===============================================================================
# modelling.sequence
#===============================================================================
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
                        Instance, Enum, List, Code, Button, Str)
from traitsui.api import (View, Item, TextEditor, CodeEditor, HGroup,
                        VGroup, Group, spring, Label, EnumEditor)
from traitsui.menu import OKButton, CancelButton

try: #hack to avoid error with py2app applications
    _scriptdir =os.path.join(os.path.dirname(__file__),'script')
    _pulse_sequence_list= os.listdir(_scriptdir)
except OSError:
    _scriptdir ='script'
    _pulse_sequence_list= os.listdir(_scriptdir)

class SequenceAvailable(HasTraits):
    type = List(_pulse_sequence_list)
sequence_available= SequenceAvailable()

class PulseSequence(HasTraits):
    
    type = Str('ideal_pulse')
    available = Constant(sequence_available)
    
    simpson_script = Code()
    script = Code()
    
    new = Button('NEW')
    dup = Button('DUPLICATE')
    rem = Button('DELETE')

    def _script_default(self):
        with open(os.path.join(_scriptdir, self.type), 'rb') as f:
            return f.read()

    def _type_changed(self):
        with open(os.path.join(_scriptdir, self.type), 'rb') as f:
            self.script = f.read()

    _name = Str()

    def _dup_fired(self):
        _dup_view = View( Item('_name', editor = TextEditor(enter_set=True),
                               label= 'name of the duplicate',),
                         buttons = [OKButton, CancelButton], width=.33)
        self._name='copy of %s'%self.type
        done = self.edit_traits(view=_dup_view, kind='modal')
        if done.result: # not cancel
            self._name = self._name.replace(' ','_')
            with open(os.path.join(_scriptdir, self._name), 'wb') as f:
                f.write(self.script)
                self.available.type.append(self._name)
            self.type = self._name

    def _new_fired(self):
        _new_view = View( Item('_name', editor = TextEditor(enter_set=True),
                               label= 'name of the new sequence',),
                         buttons = [OKButton, CancelButton], width=.33)
        self._name='new name'
        done = self.edit_traits(view=_new_view, kind='modal')
        if done.result: # not cancel
            self._name = self._name.replace(' ','_')
            with open(os.path.join(_scriptdir, self._name), 'wb') as f:
                f.write(self.script)
                self.available.type.append(self._name)
            self.type = self._name    

    def _rem_fired(self):
        pass
        #TODO

    def default_traits_view(self):
        return View(
                    Item(
                        'type',
                        editor = EnumEditor(
                                   name = 'object.available.type'),
                        label = 'sequence script'
                        ),
                    VGroup( 
                          Item(
                               'script', 
                                show_label = False
                               ),
                           ),
                    HGroup(spring,
                           Item('new', show_label=False,),
                           Item('dup', show_label=False,),
                           Item('rem', show_label=False,),
                           ),
                resizable = True,
                )

        
#===============================================================================
if __name__ == '__main__':
    PulseSequence().configure_traits()
    
    pass

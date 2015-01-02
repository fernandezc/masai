# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.core.fit.parameters
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
Model parameters handling
"""
#===============================================================================
# python import
#===============================================================================
import os
import sys
import re                       # For regular expression search
import string
import types
import logging
import importlib
import numpy as np
from UserDict import UserDict   # This is to be able to create a special dictionary

#===============================================================================
# enthought import
#===============================================================================
from traits.api import (HasTraits, Str, Instance, on_trait_change,
                        Event, Property, Button)
from traitsui.api import (View, Item, CodeEditor, HGroup)

#===============================================================================
# local import
#===============================================================================
from masai.core.fit.models import list_of_models, list_of_baselines
from masai.core.utils import DEnum

#===============================================================================
# logging
#===============================================================================
logger = logging.getLogger()

#===========================================================================
# id_generator
#===========================================================================
def _id_generator():
    """Returns a sequence of numbers for the title of the objects.
    
    Example:
    --------
    
    >>> id_generator.next()
    1
    
    """
    n=1
    while True:
        yield n
        n += 1
id_generator = _id_generator()

#===============================================================================
# FitParameters
#===============================================================================

class FitParameters(UserDict):
    """
    Allow passing a dictionary of parameters with additional properties
    to the fit function. Check if the parameter is between the specified bounds
    if any.
    
    Examples
    --------
    
    >>> fp = FitParameters()
    
    """
        
    def __init__(self):
        UserDict.__init__(self)     # Create a dictionary class
        self.lob = {}               # Lower bound
        self.upb = {}               # Upper bound
        self.fixed = {}             # true for non-variable parameter
        self.reference = {}         # 
        self.model = {}             # model to use
        self.models = []            # list of models
        self.sequence = ''          # sequence used in the experiment

    def __setitem__(self, key, value):
        key = str(key)
        if not self.reference.has_key(key):
            self.reference[key] = False
        if self.reference[key]:
            # we get a reference to another parameter 
            self.data[key] = str(value)
            self.fixed[key] = True 
        elif isinstance(value,tuple) or isinstance(value,list):
            self.data[key] = self._evaluate(value[0])
            self.lob[key] = None
            self.upb[key] = None
            try:
                if len(value)>2:
                    self.lob[key] = self._evaluate(value[1])
                    self.upb[key] = self._evaluate(value[2])
                    self._checkerror(key)
            except:
                pass
            self.fixed[key] = False
            if isinstance(value[-1],bool):
                self.fixed[key] = value[-1]
        else:
            self.data[key] = self._evaluate(value)
            self.lob[key] = None
            self.upb[key] = None
            self.fixed[key] = False

    def __getitem__(self, key):
        key = str(key)
        if self.data.has_key(key):
            return self.data[key]
        raise KeyError("parameter %s is not found"%key)

    def iteritems(self):
        return self.data.iteritems()

    def _checkerror(self, key):
        key = str(key)
        if self.lob[key] is None and self.upb[key] is None:
            return False
        elif (self.lob[key] is not None and self.data[key] < self.lob[key]) \
            or (self.upb[key] is not None and self.data[key] > self.upb[key]) :
            raise ValueError('%s value %s is out of bounds'%(key, str(self.data[key])))
        
    def __str__(self): 
        message = "#PARAMETER SCRIPT\n\nCOMMON: \n"
        for key in self.keys():
            if "_" not in  key: 
                # this is probably a common parameter 
                # (outside a model definition)
                if self.fixed[key]:
                    keystring = "\t* %s" % key
                else:
                    keystring = "\t$ %s" % key 
                lob = self.lob[key]
                upb = self.upb[key]
                if lob <= -0.1 / sys.float_info.epsilon:
                    lob = "none" 
                if upb >= +0.1 / sys.float_info.epsilon:
                    upb = "none"
                val = str(self.data[key])
                message += "%s: %s, %s, %s \n" % (keystring, val, lob, upb)
        message += ""             
        models = self.models  
        for model in models:
            message += "\nMODEL: %s\n" % model
            message += "shape: %s\n" % self.model[model]
            for key in self.keys():
                ends = "_%s" % (model)
                if not key.endswith(ends): 
                    continue 
                keystring = key[:-len(ends)] 
                
                if self.reference[key]:
                    message += "\t> %s: %s" % (keystring, self.data[key])
                    
                else:
                    if self.fixed[key]:
                        keystring = "\t* %s" % keystring
                    else:
                        keystring = "\t$ %s" % keystring    
                    lob = self.lob[key]
                    upb = self.upb[key]
                    if lob <= -0.1 / sys.float_info.epsilon:
                        lob = "none" 
                    if upb >= +0.1 / sys.float_info.epsilon:
                        upb = "none"
                    val = str(self.data[key])
                    message += "%s: %s, %s, %s" % (keystring, val, lob,  upb)
                
                message += "\n"          
        return message

    def _evaluate(self, strg):
        """
        Allow the evaluation of strings containing some operations
        
        Parameters
        ----------
        strg: string  
            A string to evaluate containing multiplier, 
            e.g., '10 k' evaluate to 10 000.
        
        Return
        ------
        value : float or bool
            Value of the string, or False, if there is an error
        
        """
        res = False
        
        if type(strg) is types.StringType or type(strg) is types.UnicodeType:
            #strg=string.upper(strg)
            p = re.compile('\s+')
            m = p.split(string.strip(strg))
            try:
                res = eval(m[0])
            except NameError:
                message="Cannot evaluate '"+strg+"' >> "+m[0]+" is not defined"
                raise NameError,message
            except SyntaxError:
                message="Syntax error in '"+strg+"'"
                raise SyntaxError,message
            #read mulitplier
            if len(m) > 1:
                try:
                    res=res*eval(m[1])
                except NameError:
                    message= "Cannot evaluate '"+strg+"' >> "+m[1]+" is not defined"
                    raise ValueError,message
                except SyntaxError:
                    message="Syntax error in '"+strg+"'"
                    raise ValueError,message
        else:
            # not a string (probably a scalar that can be return as it is)
            res = strg
            
        return res

    def to_internal(self, key):

        key = str(key)
        if not self.data.has_key(key):
            raise KeyError("parameter %s is not found"%key)

        pe = self.data[key]
        lob = self.lob[key]
        upb = self.upb[key]

        is_lob = lob is not None and lob > -0.1/sys.float_info.epsilon # lob is not None
        is_upb = lob is not None and upb < +0.1/sys.float_info.epsilon # upb is not None

        if is_lob and is_upb:
            # With min and max bounds defined
            pi = np.arcsin((2*(pe - lob)/(upb - lob)) - 1.)
        elif is_upb:
            # With only max defined
            pi = np.sqrt((upb - pe + 1.)**2 - 1.)
        elif is_lob:
            # With only min defined
            pi = np.sqrt((pe - lob + 1.)**2 - 1.)
        else:
            pi = pe
        return pi

    def to_external(self, key, pi):

        key = str(key)
        if not self.data.has_key(key):
            raise KeyError("parameter %s is not found"%key)

        lob = self.lob[key]
        upb = self.upb[key]

        is_lob = lob is not None and lob > -0.1/sys.float_info.epsilon # lob is not None
        is_upb = lob is not None and upb < +0.1/sys.float_info.epsilon # upb is not None

        if is_lob and is_upb:
            #  With min and max bounds defined
            pe = lob + ((upb - lob)/2.)*(np.sin(pi) + 1.)
        elif is_upb:
            #With only max defined
            pe = upb + 1. - np.sqrt(pi**2 + 1.)
        elif is_lob:
            #With only min defined
            pe = lob - 1. + np.sqrt(pi**2 + 1.)
        else:
            pe = pi

        self.data[key] = pe
        return pe

#===============================================================================
class FitSteps(UserDict):
#===============================================================================
    """
    Define which parameters are used in a step of the fit
    """
    count = 0
    def __init__(self, fitpar):
        UserDict.__init__(self)
        self.fixed = {}
        self.parameters = fitpar.keys

    def __setitem__(self, step, value):        
        self.data[step]=self.parameters
        self.fixed[step]={}

        # we start with all parameters fixed
        for k in self.parameters:
            self.fixed[step][k] = True

        if not value:
            value='allfixed'

        if isinstance(value,tuple) or isinstance(value,list):
            value = "|".join(value)

        if isinstance(value,str):
            if value[0]!="^": '^('+value
            if value[-1]!="$": value+")$"

            p=re.compile(value)
            for k in self.parameters:
                m = p.search(k)
                if m: self.fixed[step][k] = False
        else:
            raise ValueError, "FitStep["+str(step)+"]: Can't recognize this input!"

    def __getitem__(self, step):
        return self.fixed[step]

#===============================================================================
# ParameterScript
#===============================================================================
class ParameterScript(HasTraits):
    """
    This class allow some manipulation of the parameter list for modelling
    
    """
    from masai.plugins import IPlugin
    fd = Instance(IPlugin)
    fp = Instance(FitParameters)
        
    script = Str('')

    _list_of_models = Property(depends_on = "fd")
    modeltype = DEnum(values_name='_list_of_models')

    _list_of_baselines = Property(depends_on = "fd")
    basetype = DEnum(values_name='_list_of_baselines')

    
    #===========================================================================
    # properties
    #===========================================================================
    def _get__list_of_models(self):
        """
        Return a list of available model type, eg. voigtmodel ...
        depending on the type of experiments and the spin numbers
        """
        nd = 2 if self.fd.is_2d else 1
        spin = self.fd.get_spin()
        return list_of_models(nd, spin)

    def _get__list_of_baselines(self):
        """
        Return a list of available model type, eg. voigtmodel ...
        depending on the type of experiment
        
        """
        nd = self.fd.ndims
        return list_of_baselines(nd)

    @on_trait_change('script')
    def _check_parameters(self):
        """
        Check the validity of the parameters
        """
        self.fp = self._interpret(self.script)

    def _interpret(self, script):
        """
        Interpreter of the script content 
        """
        # init some flags
        modlabel = None
        common = False
        fixed = False
        reference = False
        
        # create a new FitParameters instance
        fp = FitParameters()
        
        # start interpreting
        lines = script.split('\n')
        lc = 0
        for item in lines:
            lc += 1  # count the lines
            line = item.strip()
            if line == '' or line.startswith("#"):
                # this is a blank or comment line, go to next line
                continue
            #split around the semi-column
            s = line.split(':')
            if len(s) != 2:
                raise SyntaxError('Cannot interpret line %d : A semi-column is missing?' % lc)
                 
            key, values = s
            key = key.strip().lower()
            if key.startswith('model'):
                modlabel = values.lower().strip()
                if modlabel not in fp.models:
                    fp.models.append(modlabel) 
                common = False
                continue
            elif key.startswith('common') or key.startswith('vars'): 
                common = True
                modlabel = 'common'
                continue
            elif key.startswith('shape'):
                shape = values.lower().strip()
                if shape is None or (shape not in self._list_of_models and shape not in self._list_of_baselines):
                    raise SyntaxError('Shape of this model "%s" was not specified or is not implemented'%shape)
                fp.model[modlabel] = shape
                common = False
                continue
            else:
                if modlabel is None and not common:
                    raise SyntaxError('The first definition should be a label for a model or a block of variables or constants.')                      
                # get the parameters
                if key.startswith('*'):
                    fixed = True
                    reference = False
                    key = key[1:].strip() 
                elif key.startswith('$'):
                    fixed = False
                    reference = False
                    key = key[1:].strip()
                elif key.startswith('>'):
                    fixed = True
                    reference = True
                    key = key[1:].strip()
                else:
                    raise SyntaxError('Cannot interpret line %d: A parameter definition must start with *,$ or >' % lc)
                    
                # store this parameter
                s = values.split(',')
                if len(s) > 3:
                    raise SyntaxError('line %d: value, min, max shoud be defined in this order' % lc)
                    s = s[:3]
                elif len(s) == 2:
                    raise SyntaxError('only two items in line %d' % lc)
                    s.append('none')
                elif len(s) == 1:
                    s.extend(['none', 'none'])
                value, mini, maxi = s
                if mini.strip().lower() in ['none', '']: mini = str(-1. / sys.float_info.epsilon)
                if maxi.strip().lower() in ['none', '']: maxi = str(+1. / sys.float_info.epsilon)
                if modlabel != 'common':
                    ks = "%s_%s" % (key, modlabel)
                else:
                    ks = "%s" % (key)  
                fp.reference[ks] = reference
                if not reference:
                    fp[ks] = value.strip(), mini.strip(), maxi.strip(), fixed      
                else:
                    fp[ks] = value.strip()
        return fp 

                
if __name__ == '__main__':

    fp = FitParameters()
    fp['one']  = 1, -1, 3, True

    pass
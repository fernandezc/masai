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
from .models import list_of_models, list_of_baselines
from masai.core.utils import DEnum

#===============================================================================
# logging
#===============================================================================
logger = logging.getLogger()

#===============================================================================
# constants
#===============================================================================
FIXED=True

doc_params = """
#-----------------------------------------------------------
# syntax for parameters definition:
# name: value, low_bound,  high_bound
#  * for fixed parameters
#  $ for variable parameters
#  > for reference to a parameter in the COMMON block
#    (> is forbidden in the COMMON block)
# common block parameters should not have a _ in their names
#-----------------------------------------------------------
"""

default_params_script = """
COMMON:
# common parameters ex.
# $ gwidth: 1.0, 0.0, none
# $ gratio: 0.5, 0.0, 1.0

"""
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
#TODO: Traitify
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
            raise ValueError, "FitStep["+str(step)+"]: Can't reconize this input!"

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
        
    script_filename = Str('')
    script = Str('')

    _list_of_models = Property(depends_on = "fd")
    modeltype = DEnum(values_name='_list_of_models')
    addmodel = Button('ADD LINE')

    _list_of_baselines = Property(depends_on = "fd")
    basetype = DEnum(values_name='_list_of_baselines')
    addbase = Button('ADD BASELINE')

    updated = Event   
    
    lastparfailure = Instance(Exception,allow_none=True) 

    
    #===========================================================================
    # properties
    #===========================================================================
    #@cached_property
    def _get__list_of_models(self):
        """
        Return a list of available model type, eg. voigtmodel ...
        depending on the type of experiments and the spin numbers
        """
        nd = self.fd.ndims
        spin = self.fd.NUC1[-1].spin
        return list_of_models(nd, spin)

    #@cached_property
    def _get__list_of_baselines(self):
        """
        Return a list of available model type, eg. voigtmodel ...
        depending on the type of experiment
        
        """
        nd = self.fd.ndims
        return list_of_baselines(nd)
    
    #===========================================================================
    # view
    #===========================================================================
    def default_traits_view(self):        
        return View(
#                Group( 
                      Item('content',
                                  show_label=False,
                                  style='custom',
                                  #height=400,
                                  editor=CodeEditor(auto_set=True),
                        ),
                        HGroup(
                              Item('addmodel', 
                                   show_label=False,),
                              Item('modeltype', 
                                   label='model type',),
                              visible_when='modeltype is not None and not lastparfailure'
                        ),
                        HGroup(
                              Item('addbase', 
                                   show_label=False,),
                              Item('basetype', 
                                   label='baseline',),
                              visible_when='basetype is not None and not lastparfailure'
                        ),
                        Item('lastparfailure',
                             show_label=False, 
                             visible_when='lastparfailure', 
                             style='readonly', 
                             emphasized=True),
#                    ),
                    resizable=True,
                )
        
    @on_trait_change('fd')
    def _fdchanged(self):
        """
        This event is fired when the source is first defined or changed
        """
        logger.debug("initialisation of parameter class")
        
        #try to load a previously open parameter script file
        self.scriptfilename = os.path.join(self.fd.path, "param.fit")
        if os.path.exists(self.scriptfilename):
            #if it is found, set content  
            with open(self.scriptfilename, 'rb') as f:
                self.content = f.read()
        else:
            #else, we load some default script
            self.content = default_params_script
        
        # define a timer
        if not self.timer:
            self.timer = Timer(3.0, self._check_parameters)
            self.timer.start()
            
    def _check_parameters(self):
        """
        Check the validity of the parameters
        """
        try:
            self.lastparfailure = None
            if self.content != self.savecontent:
                logger.debug("parameters: content updated")
                self.fp = self._interpret(self.content)
                self._save_parameters()
                self.updated = True
            self.savecontent = self.content
            return True
        except Exception,e:
            self.lastparfailure = e
            logger.error(e)
            return False
        
    def _save_parameters(self):
        """
        Save the content to a file
        """   
        if self.scriptfilename and not self.lastparfailure:  
            logger.debug("parameters: content saved")
            f = open(self.scriptfilename, 'wb')
            f.write(self.content)
            f.close()            

#    @on_trait_change('content')
#    def _update_parameters(self, name, old, new):
#        """
#        fired if some changes are detected in the script content
#        """        
#        if self._check_parameters():
#            logger.debug("parameters: content updated")
#            self._save_parameters()   
         
    def _addmodel_fired(self, baseline=False):
        """
        Add action fired. We import the required model and some parameters
        before writing in the script content
        """
        self.fd.show_model=True
        
        module = importlib.import_module("masai.modelling.models")
        if not baseline:
            class_ = getattr(module, self.modeltype)
        else:
            class_ = getattr(module, self.basetype)
        script = class_().script
        
        # we determine some parameters from the source data
        uc = self.fd.units
        x = uc[-1].ppm_scale() # scale in ppm
        y = self.fd.data.real # np.abs(self.fd.data[lx/4:-lx/4])
        try:
            ym = self.fd.modeldata[-1]
            y = y-ym
        except:
            pass
        ampl = np.max(y)
        argm = np.argmax(y) #+lx/4
        pos = x[argm]
        poslb, poshb = uc[-1].ppm_limits()
        w= np.argwhere(y > ampl/2.)
        wmin = np.min(w) #+lx/4
        wmax = np.max(w) #+lx/4
        xm = min( abs(x[wmax] - pos), abs(x[wmin] - pos))
        width = xm * 2. # sigma corresponding to the FWHM
#        #correct the ampl to take into account the width
        delta = abs(x[1]-x[0])
        ampl = ampl * np.pi * width /2. /delta
        newid = id_generator.next()
        #try to get the amplitude for polynomial baseline
        a = y[0]
        b = y[np.argmin(y)]
        scale = abs(a-b)/abs(x[0]-x[np.argmin(y)])**2/ampl/10.
        name = 'line'
        if baseline:
            name='baseline'
        while '%s%d'%(name,newid) in self.fp.models:
            newid = id_generator.next()
            
        #norm the amplitude
        ampl = ampl/self.fd.norm
        self.content += "\n" + script % {"id":newid, "width":width, "ampl":ampl, 
                         "pos":pos, "poslb":poslb, "poshb":poshb, 'scale':scale}

    def _addbase_fired(self):
        """ Add baseline action fired. We import the required model and some parameters
        before writing in the script content
        """
        self._addmodel_fired(baseline=True)
        
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

                

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
#==============
# python import
#==============
import sys
import re                       # For regular expression search
import string
import types
import numpy as np
from UserDict import UserDict # This is to be able to create a special dictionary

#=================
# enthought import
#=================
from traits.api import (HasTraits, Str, Instance, on_trait_change,
                        Property, Int, Float, List)

#=============
# local import
#=============
from masai.core.fit.models import list_of_models, list_of_baselines
from masai.core.utils import DEnum
#from masai.core.source import Source

#=============
# id_generator
#=============
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


#==============
# FitParameters
#==============
class FitParameters(UserDict):
    """
    Allow passing a dictionary of parameters with additional properties
    to the fit function. Check if the parameter is between the specified bounds
    if any.
    """

    #------------------
    def __init__(self):
        UserDict.__init__(self)     # Create a dictionary class
        self.lob = {}               # Lower bound
        self.upb = {}               # Upper bound
        self.fixed = {}             # true for non-variable parameter
        self.reference = {}         #
        self.common = {}            # indicates if a parameters belong to a common block
        self.model = {}             # model to use
        self.models = []            # list of models
        self.sequence = ''          # sequence used in the experiment
        self.expvars = []           # list of parameters which are experiment dependent
        self.expnumber = 1          # number of experiments

    #---------------------------------
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

    #--------------------------
    def __getitem__(self, key):
        key = str(key)
        if self.data.has_key(key):
            return self.data[key]
        raise KeyError("parameter %s is not found"%key)

    #-------------------
    def iteritems(self):
        return self.data.iteritems()

    #--------------------------
    def _checkerror(self, key):
        key = str(key)
        if self.lob[key] is None and self.upb[key] is None:
            return False
        elif (self.lob[key] is not None and self.data[key] < self.lob[key]) \
            or (self.upb[key] is not None and self.data[key] > self.upb[key]) :
            raise ValueError('%s value %s is out of bounds'%(key, str(self.data[key])))

    #-----------------
    def __str__(self): 

        #.............................................................
        def makestr(key):

            keystring = key.split('_')[0]
            if self.reference[key]:
                return "\t> %s: %s \n" % (keystring, self.data[key])
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

                return "%s: %s, %s, %s \n" % (keystring, val, lob, upb)
        #..............................................................

        message = "#PARAMETER SCRIPT\n\nCOMMON: \n"

        var = ""
        for item in self.expvars:
            var += " %s"%item

        if var:
            message += "\texperiment_number: %s \n"%str(self.expnumber)
            message += "\texperiment_variables: %s \n"%var

        # look for common parameters
        for key in self.keys():
            keysp = key.split('_')[0]
            if self.common[keysp]:
                message += makestr(key)

        # model parameters
        models = self.models  
        for model in models:
            message += "\nMODEL: %s\n" % model
            message += "shape: %s\n" % self.model[model]
            for key in self.keys():
                keyspl = key.split('_')
                if not model in '_'.join(keyspl[1:]):
                    continue
                message += makestr(key)
        return message

    #-------------------------
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

    #----------------------------------------
    def to_internal(self, key, expi=None):
        """
        if expi is not none, several parameters to create
        """
        key = str(key)
        if not self.data.has_key(key):
            raise KeyError("parameter %s is not found"%key)

        if expi is not None:
            pe = self.data[key][expi]
        else:
            pe = self.data[key]
        lob = self.lob[key]
        upb = self.upb[key]

        is_lob = lob is not None and lob > -0.1/sys.float_info.epsilon # lob is not None
        is_upb = lob is not None and upb < +0.1/sys.float_info.epsilon # upb is not None

        if is_lob and is_upb:
            lob = min(pe,lob)
            upb =max(pe,upb)
            # With min and max bounds defined
            pi = np.arcsin((2*(pe - lob)/(upb - lob)) - 1.)
        elif is_upb:
            upb =max(pe,upb)
            # With only max defined
            pi = np.sqrt((upb - pe + 1.)**2 - 1.)
        elif is_lob:
            lob = min(pe,lob)
            # With only min defined
            pi = np.sqrt((pe - lob + 1.)**2 - 1.)
        else:
            pi = pe
        return pi

    #------------------------------
    def to_external(self, key, pi):

        key = str(key)
        if not self.data.has_key(key):
            raise KeyError("parameter %s is not found"%key)

        lob = self.lob[key]
        upb = self.upb[key]

        is_lob = lob is not None and lob > -0.1/sys.float_info.epsilon # lob is not None
        is_upb = lob is not None and upb < +0.1/sys.float_info.epsilon # upb is not None

        if not isinstance(pi, list):
            pi = [pi, ]  #make a list

        pe = []
        for item in pi:
            if is_lob and is_upb:
                #  With min and max bounds defined
                pei = lob + ((upb - lob)/2.)*(np.sin(item) + 1.)
            elif is_upb:
                #With only max defined
                pei = upb + 1. - np.sqrt(item**2 + 1.)
            elif is_lob:
                #With only min defined
                pei = lob - 1. + np.sqrt(item**2 + 1.)
            else:
                pei = pi
            pe.append(pei)

        if len(pe)==1:
            pe = pe[0]

        self.data[key] = pe

        return pe

#========================
class FitSteps(UserDict):
#========================
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

#================
# ParameterScript
#================
class ParameterScript(HasTraits):
    """
    This class allow some manipulation of the parameter list for modelling
    
    """
    from masai.plugins import IPlugin
    sources = List(Instance("masai.plugins.sources.bruker_plugin.Bruker"))
    fp = Instance(FitParameters)
        
    script = Str('')

    _list_of_models = Property()
    modeltype = DEnum(values_name='_list_of_models')

    _list_of_baselines = Property()
    basetype = DEnum(values_name='_list_of_baselines')

    
    #===========================================================================
    # properties
    #===========================================================================

    #------------------------------
    def _get__list_of_models(self):
        """
        Return a list of available model type, eg. voigtmodel ...
        depending on the type of experiments and the spin numbers
        """
        models = []
        for source in self.sources:
            spin = source.get_spin()
            n_dims = source.ndims
            models.extend(list_of_models(n_dims, spin))
        return set(models)

    #---------------------------------
    def _get__list_of_baselines(self):
        """
        """
        baselines = []
        for source in self.sources:
            n_dims = source.ndims
            baselines.extend(list_of_baselines(n_dims))
        return set(baselines)

    #--------------------------
    @on_trait_change('script')
    def _check_parameters(self):
        """
        Check the validity of the parameters
        """
        self.fp = self._interpret(self.script)

    #-----------------------------
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

        # set the number of experiments
        fp.expnumber = len(self.sources)
        print("The number of experiment is set to %d"%fp.expnumber)

        # start interpreting ------------------------------------------------------
        lines = script.split('\n')
        lc = 0

        for item in lines:
            lc += 1  # -------------- count the lines
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
            elif key.startswith("experiment"): # must be in common
                if not common:
                    raise SyntaxError(
                    "'experiment_...' specification was found outside the common block.")
                if "variables" in key:
                    expvars = values.lower().strip()
                    expvars = expvars.replace(',',' ').replace(';',' ')
                    expvars = expvars.split()
                    fp.expvars.extend(expvars)
                continue
            else:
                if modlabel is None and not common:
                    raise SyntaxError(
                   "The first definition should be a label for a model or a block of variables or constants.")
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
                s = map(string.strip, s)
                if len(s)>1 and ('[' in s[0]) and (']' in s[1]): #list
                    s[0] = "%s, %s"%(s[0], s[1])
                    if len(s)>2:
                        s[1:] = s[2:]
                if len(s) > 3:
                    raise SyntaxError('line %d: value, min, max should be defined in this order' % lc)
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
                    fp.common[key] = False
                else:
                    ks = "%s" % (key)
                    fp.common[key] = True
                #if key in fp.expvars:
                #    for i in xrange(len(self.sources)):
                #        ks = "%s_exp%d"%(ks, i)
                fp.reference[ks] = reference
                if not reference:
                    val = value.strip()
                    val = eval(val)
                    if isinstance(val, list):
                        # if the parameter is already a list, that's ok if the number of parameters is ok
                        if len(val)!=fp.expnumber:
                            raise ValueError('the number of parameters for %s is not the number of experiments.'%len(val))
                        if key not in fp.expvars:
                            raise ValueError('parameter %s is not declared as variable'%key)
                    else:
                        if key in fp.expvars:
                            # we create a list of parameters corresponding
                            val = [val]*fp.expnumber
                    fp[ks] = val, mini.strip(), maxi.strip(), fixed
                else:
                    fp[ks] = value.strip()


        return fp 

###########################

if __name__ == '__main__':

    import os, commands
    from masai.api import Bruker, Process, Plot1d, Fit1d

    # we try here a multi source fit
    # print(commands.getoutput('pwd'))

    DATADIR = os.path.expanduser(os.path.join('..', '..', '..', 'exemples'))
    user = 'user'
    name = 'RUB-17'

    path1 = os.path.join(DATADIR, user, 'nmr', name, str(1))
    source1 = Bruker(path1)
    Process(source1, transform='em 1; zf 0; ft')

    path2 = os.path.join(DATADIR, user, 'nmr', name, str(2))
    source2 = Bruker(path2)
    Process(source2, transform='em 1; zf 0; ft')

    #p = Plot1d(source1, hold=True)
    #p = Plot1d(source2, axe=p.axe)

    # fit
    # try with variable width
    pos = [-94.9200018156, -95.0999986491, -94.3173715978, -93.9072128893,
           -92.8999998621, -92.6363220634, -82.0826836442]

    script = """
    COMMON:
      experiment_variables: ampl, width
      $ gratio: 0.5, 0.0, 1.0

    """
    for i, delta in enumerate(pos):
        fa = '*' if i==0 else '$'
        script +="""
            MODEL: LINE_%d
            shape: voigtmodel
                %s ampl:  0.6, 0.0, none
                $ pos: %.3f, %.3f, %.3f
                > ratio: gratio
                $ width: [0.4, %.3f], %.3f, %.3f
        """%(i, fa, delta, delta-.2, delta+0.2,
                    0.57, 0.45, 0.65)

    print('#'*10)
    # prepare fit
    f = Fit1d( [source1, source2], script)
    f.run(maxiter=200, every=10)
    print('#'*10)
    print(f.script)
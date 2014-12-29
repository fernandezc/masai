#
#===============================================================================
# controllers.fitcontroller
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
"""
.. _fitcontroller:

This module defines objets and methods related to the modelling 
and fitting procedure of the datasets

 """
from traits.api import (HasTraits, Str, Instance, on_trait_change,
                        Enum, Float, Int, Callable, Property,
                        cached_property, Button, Bool, Array, List,
                        DelegatesTo)
from traitsui.api import (View, Item, Label, TextEditor,
                          Group, HGroup, VGroup, spring,)
import os.path
import sys
import numpy as np
import re
from threading import Thread

from masai.core.fit.parameters import ParameterScript
from masai.core.fit.models import getmodel
from masai.core.fit.sequence import PulseSequence
#from masai.api import findregion
from masai.core.fit.optimize import optimize # a modified copy of scipy.optimize
         
import logging
logger = logging.getLogger()

doc_params = """
#-----------------------------------------------------------
# syntax for parameters definition:
# name: value, low_bound,  high_bound
# *for fixed parameters
#  $ for variable parameters
#  > for reference to a parameter in the COMMON block
#    (> is forbidden in the COMMON block)
# common block parameters should not have a _ in their names
#-----------------------------------------------------------
"""

default_params_script = """
COMMON:
$ globgb: %(gb).3f, 0.0, none
$ globlb: 0.001, 0.0, none
$ globamp: %(ampl).3f, 0.0, none

MODEL: Line1
shape: voigtmodel
* amp:  1.0, 0.01, none
$ pos:  %(pos).3f, %(poslb).3f, %(poshb).3f
> gb:  globgb
> lb:  globlb

"""

doc_model = """
#-----------------------------
# TODO: make the documentation
#-----------------------------
"""

head_model = """
def userfunction(modeldata, pars):
 
#*******************************
"""

default_userscript = ""

tail_model = \
"""
#*******************************
    return modeldata
"""

#===============================================================================
# global function for modelfunctions
#===============================================================================

def parsing(expr, param):
    keyword = r"\b([a-z]+)\b" #match a whole word
    pat = re.compile(keyword)
    mo = pat.findall(expr)
    for kw in mo:
        if param.has_key(kw):
            expr = expr.replace(kw, str(param[kw]))
        elif np.__dict__.has_key(kw):# check if it is a recognized math function
            expr = expr.replace(kw, "np.%s" % kw)
    return expr
    
def prepare(param):
    """
    This function is needed for the script related to modelfunction
    """
    new_param = param.copy()
    
    for key in param.iterkeys():
        if param.reference[key]:
            new_param.reference[key] = False    # important to put it here 
                                                # before other instruction
            # try to interpret the given refpar expression
            refpar = param[key]
            while True:
                par = parsing(refpar, new_param)
                if par == refpar: break
                refpar = par

            try:                                    
                new_param[key] = eval(str(refpar))
            except:
                raise ValueError('Cannot evaluate the expression: %s: %s'
                                 % (key, param[refpar]))
                return
            new_param.fixed[key] = True
            new_param.reference[key] = True     #restore it for the next call 

    return new_param

#===============================================================================
# makefunc
#===============================================================================
def _makefunc(script=default_userscript):
    
    codescript = head_model
    scriptlines = script.split('\n')
    for line in scriptlines:
        codescript += "\t%s\n" % line
    codescript += tail_model
    
    #logger.debug(codescript)
    
    userfunction = None
    #try:
    code = compile(codescript, '<script>', 'exec')
    exec(code)
    return userfunction
    
    #except:
    #    return None

#===============================================================================
# automatic calculation of amplitude and baseline
#===============================================================================
def ampandbas(xi, expe, calc):
    """Automatically calculate correct amplitude A and baseline 
    (baseline linear model a*i+b) by detemining the zero of the first derivative
     with respect to A, a, and b
    """
    n = xi.size
    sic = sum(xi*calc)
    sc = sum(calc)
    sie = sum(xi*expe)
    se = sum(expe)
    sce = sum(expe*calc)
    sc2 = sum(calc*calc)
    A = (
         (12*n*sic + (-6*n**2 + 18*n - 12)*sc)*sie
         + (-6*n**2 + 18*n - 12)*se*sic 
         + (4*n**3 - 6*n**2 + 2*n)*sc*se
         + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sce)/(12*n*sic**2
                                                + (-12*n**2 + 36*n - 24)*sc*sic
         + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2 
         + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    a = -(
          (12*n*sc2 - 12*sc**2)*sie + (12*sc*se - 12*n*sce)*sic 
          + (-6*n**2 + 18*n - 12)*sc2*se
          + (6*n**2 - 18*n + 12)*sc*sce)/(12*n*sic**2 
                                                + (-12*n**2 + 36*n - 24)*sc*sic
          + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2 
          + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    b = -(
          (12*sc*sic + (-6*n**2 + 18*n - 12)*sc2)*sie
          - 12*se*sic**2 + (6*n**2 - 18*n + 12)*sce*sic 
          + (4*n**3 - 6*n**2 + 2*n)*sc2*se
          + (-4*n**3 + 6*n**2 - 2*n)*sc*sce)/(12*n*sic**2
          + (-12*n**2 + 36*n - 24)*sc*sic 
          + (-n**4 - 12*n**3 + 37*n**2 - 36*n + 12)*sc2
          + (4*n**3 - 6*n**2 + 2*n)*sc**2)
    # in case the modeldata is zero, to avoid further errors
    if np.isnan(A): A=0.0
    if np.isnan(a): a=0.0
    if np.isnan(b): b=0.0
    return A, a, b    

#===============================================================================
# ModelCompute
#===============================================================================
class ModelComputeThread(Thread):
    """Inspired from gael varoquaux tutorial
    
    """
    def run(self):
        control = self.control
        if not control.sequence:
            return  # not fully initialized yet, return
        try:
            logger.debug(str(self)+' started')
            control.lastfitfailure = None
            modeldata, names = control._get_modeldata() 
        except Exception, e:
            control.lastfitfailure = e
            control.paramscript.lastparfailure = e
            logger.error(e)
            logger.debug(str(self)+' stopped on error')
            return
            
        self.control.fd.modelnames = names   # must come first 
                                            #(to avoid plotting routine errors)
        self.control.fd.modeldata = modeldata
        logger.debug(str(self)+' finished')
        
    def stop(self):
        logger.debug(str(self)+' stopped')
        self._Thread__stop()

#===============================================================================
# ModelFitting
#===============================================================================
class ModelFittingThread(Thread):
    """Inspired from gael varoquaux tutorial
    """
    stopped = False
    def run(self):
        control = self.control
        try:
            control.lastfitfailure = None
            control._fitting()
        except Exception, e:
            #some common errors
            if e.message.strip() == \
                       "zero-size array to maximum.reduce without identity":
                control.lastfitfailure = ValueError(
                           "At least one parameter should be fitted\n"\
                           + "(use $ in front of the parameter to fit)")
            else:
                control.lastfitfailure = e
            logger.error(e)
            
        control._running = False
             
#===============================================================================
# Class FitController
#===============================================================================
class FitController(HasTraits):
    
    fd = Instance('masai.fileio.source.Source')
    
    paramscript = Instance(ParameterScript) 
    userscriptfile = Str('')
    userscriptcontent = Str
    userfunction = Property(Callable,
                            depends_on='userscriptcontent, fd.is_resized')
    lastfitfailure = Instance(Exception, allow_none=True)
    computethread = Instance(ModelComputeThread)
    fittingthread = Instance(ModelFittingThread)
    
    sequence = Instance(PulseSequence)
    
    startbtn = Button('START FIT')
    stopbtn = Button('STOP FIT')

    fp = Property(depends_on='paramscript.updated')
    fittable = DelegatesTo('fd')
    
    method = Enum(("SIMPLEX"))#, "BFGS", "CG"))
    tolerance = Float(1.e-5)
    maxiter = Int(250)
    maxfuncalls = Int(2500)
    plotevery = Int(10)
    
    funcalls = Int
    iterations = Int
    fopt = Float
    warnmess = Str
    fitpar = Str
    
    intensities = Str('')
    A = DelegatesTo('fd')
        
    #===========================================================================
    # private
    #===========================================================================        
    _loading = Bool(True)
    _running = Bool(False)
    
    _save_modeldata = Array(allows_none=True)
    _save_names= List
    _save_fp = Instance('masai.modelling.parameters.FitParameters')     
      
    #===========================================================================
    # default_traits_view
    #===========================================================================
    def default_traits_view(self):        
        return View(
                VGroup(
                    Group(
                        Group(Item('paramscript',
                                    show_label=False,
                                    style='custom',
                                    #enabled_when='not _running'
                                ),
                               label="Parameters",
                               dock='tab',
                        ),
                        Group(
                              Item(name='intensities',
                                   show_label = False,
                                   style = 'readonly',
                                   editor=TextEditor(multi_line=True)
                                   ),
                              label= 'Report',
                              dock = 'tab',
                              ),
                        #=======================================================
                        ## NOT USED FOR NOW
                        # Group( Item('userscriptcontent',
                        #             show_label=False,
                        #             style='custom',
                        #             editor=CodeEditor(lexer='python', 
                        #                               auto_set=True)),
                        #       label="User script",
                        #       dock='tab',
                        # ),
                        # Group( Item('sequence',
                        #             style='custom',
                        #             show_label =False),
                        #       label="pulse sequence",
                        #       dock='tab',
                        # ),
                        #=======================================================
                        springy=True,
                        layout='tabbed',
                    ),
                    Group(
                        VGroup(
                               #Item('method', label='Method',), #not used for now
                               Item('tolerance', width=50),
                               Item('maxiter',
                                    label="Max of iterations",
                                    width=50),
                               Item('maxfuncalls',
                                    label='Max of calls/iterations',
                                    width=50),
#                               Item('plotevery',
#                                    label='Nb of iterations/plot updates',
#                                    width=50),
                               orientation='vertical',
                        ),
                        VGroup(
                               HGroup(
                               Item('iterations',
                                    label=u"iterations",
                                    style='readonly',
                                    width=50),
                               Item('funcalls',
                                    label='calls',
                                    style='readonly',
                                    width=50),
                               Item('fopt',
                                    label='Error (Chi)',
                                    style='readonly',
                                    width=100),
                               ),
                               Item('warnmess',
                                    editor=TextEditor(multi_line=True),
                                    show_label=False,
                                    width=300,
                                    style='readonly',
                                    emphasized=True,
                                    visible_when="funcalls>0",),
                               Item('lastfitfailure',
                                    editor=TextEditor(multi_line=True),
                                    show_label=False,
                                    visible_when='lastfitfailure',
                                    style='readonly',
                                    width=350,
                                    emphasized=True),
                               orientation='vertical',
                               show_border=True,
                               label='results',
                        ),
                        Item('_'),
                        HGroup(
                              Item('startbtn', show_label=False,
                                    visible_when='not _running',
                                    ),
                               Item('stopbtn', show_label=False,
                                    visible_when='_running'),
                               spring,
                               orientation='horizontal',
                        ),
                        orientation='vertical',
                    ),
                    visible_when='fittable',
                ),
                 VGroup(
                   Item('_'),
                   Label('WARNING'),
                   Label('A fourier transform must be done to access fitting'),
                   Item('_'),
                   visible_when='not fittable or fd is None',
                ),
                resizable=True,
                scrollable=True,
            )
        
    #===========================================================================
    # fd event router
    #===========================================================================
    @on_trait_change('fd.+')
    def _eventrouter(self, objet=None, name=None, old=None, new=None):
        
        if not self.fd.is_loaded: 
            # we return if the source is not valid
            return

        if name == 'fd':
            logger.debug("initialisation of Fit class")

#            self.model = ''
#            self.simpson = ''

            self.paramscript = ParameterScript(fd=self.fd)

#            # set userscriptcontent       
            if not os.path.exists(self.userscriptfile):
                self.userscriptcontent = default_userscript
            else:
                self.userscriptfile = os.path.join(self.fd._directory, 
                                               "userscript.fit")
                with open(self.userscriptfile, 'rb') as f:
                    self.userscriptcontent = f.read()
        
            _type = self.fd.sourcepref.get("fit.sequence", 'ideal_pulse')  
            self.sequence = PulseSequence(type=_type)


    @cached_property
    def _get_fp(self):
        return self.paramscript.fp
           
    def _save_userscript(self):
        if self.userscript:
            f = open(self.userscriptfile, 'wb')
            f.write(self.userscriptcontent)
            f.close()           
    
    @on_trait_change('userscriptcontent')
    def _update_userscriptcontent(self, name, old, new):
        self._save_userscript()
        self.updated = True
          
    @cached_property
    def _get_userfunction(self): 
        """user function 
        """
        #logger.debug('get userfunction')
        #create or update the user function
        userfunction = _makefunc(self.userscriptcontent)
        return userfunction

    #===========================================================================
    # _get_modeldata
    #===========================================================================
    def _get_modeldata(self):
        """Calculate the model
        """
        
        if not self.fd.forceupdatemodel and \
                         hasattr(self, '_save_fp') and self.fp == self._save_fp:
            # absolutely no changes, just return the saved modeldata
            return self._save_modeldata, self._save_names
        
        self.fd.forceupdatemodel=False # suppress force update
        self._save_fp = self.fp.copy()

        # Prepare parameters
        parameters = prepare(self.fp)

        # Get the list of models
        models = self.fp.models
        nbmodels = len(models)

        # Make an array 'modeldata' with the size of the source of data 
        # which will contains the data produced by the models
        # This name must always be 'modeldata' 
        # which will be returned to the main program.
        
        uc = self.fd.units
        x = uc[-1].ppm_scale()
        modeldata = np.zeros((nbmodels + 2, x.size)).astype(float)
    
        if nbmodels < 1:
            names = ['baseline', 'modelsum']
            self._save_modeldata = modeldata
            self._save_names = names
            return modeldata, names

        # Calculates model data
        # The first row (i=0) of the modeldata array is the baseline, 
        # so we fill the array starting at row 1
        row = 0
        names = ['baseline', ]
        for model in models:
            if not model.startswith('baseline'):
                row += 1
                modeldata[row] = getmodel(x, model, parameters,
                                          fd=self.fd, sequence=self.sequence)
                names.append(model)
            else:
                modeldata[0] += getmodel(x, model, parameters,
                                         fd=self.fd, sequence=self.sequence)
                
        # call the userfunction 
        # (in case something is executed in it - nothing by default)
        if default_userscript:
            modeldata = self.userfunction(modeldata, self.fp)  
            
        # make the sum
        modeldata[row + 1] = modeldata.sum(axis=0)
        names.append('modelsum')

        #remove unused row
        modeldata = modeldata[:row + 2]
        
        #calculate automatic baseline and amplitude
        data = self.fd.data.real
        xi = np.arange(x.size)
        
        # take into account fitzone
        uc = self.fd.units
        x = uc[-1].ppm_scale()
        fitzone = findregion(x, self.fd.zone)

        A, a, b = ampandbas(fitzone-fitzone[0], data.take(fitzone), 
                            modeldata[-1].take(fitzone))
        
        modeldata = A*modeldata
        baseline = a*(xi-fitzone[0]) + b
        
        #update the modeldata
        modeldata[0] += baseline
        modeldata[-1] += baseline
        
        self._save_modeldata = modeldata
        self._save_names = names
        self.A = A
        
        # calculate intensities
        at = 0.
        for model in models:
            if model.startswith('baseline'):
                continue
            at += parameters["ampl_%s"%model]
            
        self.intensities = 'Intensities:\n       \tppm   \tRelative (%)     \tAbsolute\n'
        for model in models:
            if model.startswith('baseline'):
                continue
            self.intensities += \
                '%(model)7s:\t%(pos).2f\t%(ampl)g \t%(conc)g\n'%{
                                  'model':model,
                                  'pos':parameters["pos_%s"%model],
                                  'ampl':parameters["ampl_%s"%model]*100./at,
                         'conc':parameters["ampl_%s"%model]*self.A/self.fd.norm}

        return modeldata, names
    
    #===========================================================================
    # recompute models
    #===========================================================================
    @on_trait_change('fp, fd.is_resized, fd.forceupdatemodel')
    def _updated_fp(self, objet=None, name=None, old=None, new=None):
                
        self.lastfitfailure = None
        if hasattr(self.paramscript,'lastparfailure') \
                   and self.paramscript.lastparfailure: 
            return
        
        #check is the thread is already running. If yes send abort first
        if self.computethread and self.computethread.isAlive():
            self.computethread.stop()
        self.computethread = ModelComputeThread()
        self.computethread.control = self
        self.computethread.start()
        
        self.fd.model_need_replot

    #===========================================================================
    # Fitting actions
    #===========================================================================
    def _startbtn_fired(self):
        """
        Start fitting 
        """
        logger.info('\nfit started...')
        logger.info(
                'Take a cup of coffee and be patient as it may take a while!')
        
        self.fd.can_plot=False
        self.funcalls = 0
        self._running = True
        self.warnmess = ''
        #self.lastfitfailure = ''
        
        if len(self.fd.zone) < 1:
            self._running = False
            return  # this should not happen
 
        
        #check is the thread is already running. If yes send abort first
        if self.fittingthread and self.fittingthread.isAlive():
            self.fittingthread.wants_abort = True
        else:
            self.fittingthread = ModelFittingThread()
            self.fittingthread.wants_abort = False
            self.fittingthread.control = self
            self.fittingthread.start()       

    def _stopbtn_fired(self):
        """
        Stop fit  : TODO make it working using Threads
        """
        if self.fittingthread and self.fittingthread.isAlive():
            self.fittingthread.stopped = True
        self.fd.can_plot=True
        #logger.info('fit stopped. Wait for the last execution call')
        return True
                
    #===========================================================================
    # _fitting
    #===========================================================================

    #  Internal/external transformation
    #  These transfomations are used in the MINUIT package, 
    #  and described in detail 
    #  in the section 1.3.1 of the MINUIT User's Guide.
    
    def internal(self, pe, lob, upb):
        islob = (lob > -0.1/sys.float_info.epsilon) 
        isupb = (upb < +0.1/sys.float_info.epsilon)
        
        if islob and isupb:
            # With min and max bounds defined
            pi = np.arcsin((2*(pe - lob)/(upb - lob)) - 1.)
        elif isupb:
            # With only max defined
            pi = np.sqrt((upb - pe + 1.)**2 - 1.)
        elif islob:
            # With only min defined
            pi = np.sqrt((pe - lob + 1.)**2 - 1.)
        else:
            pi = pe    
        return pi
    
    def external(self, pi, lob, upb):
        islob = (lob > -0.1/sys.float_info.epsilon) 
        isupb = (upb < +0.1/sys.float_info.epsilon)

        if islob and isupb:
            #  With min and max bounds defined
            pe = lob + ((upb - lob)/2.)*(np.sin(pi) + 1.)
        elif isupb:
            #With only max defined
            pe = upb + 1. - np.sqrt(pi**2 + 1.)
        elif islob:
            #With only min defined
            pe = lob - 1. + np.sqrt(pi**2 + 1.)
        else:
            pe = pi
        return pe
 
    def _fitting(self):   
        """Main fitting procedure
        """ 
        logger.info('Entering fitting procedure') 
        self._forcereplot = False
        
        uc = self.fd.units
        x = uc[-1].ppm_scale()
        fitzone = findregion(x, self.fd.zone)
                 
        # internally defined function chi2
        def chi2(params, keys):
            """ 
            Return sum((y - x)**2) 
            
            """
            for i in range(len(params)):
                key = keys[i]
                # restore external parameters in case of bounding
                pe = self.external(params[i],
                                   self.fp.lob[key], self.fp.upb[key])
                self.fp[key] = pe, self.fp.lob[key], self.fp.upb[key]

            # model spectrum
            modeldata = self._get_modeldata()[0]
                
            # baseline is already summed with modeldata[-1]
            diff = (self.fd.data.real - modeldata[-1])**2 
            chi2 = np.sum(diff.take(fitzone))

            return chi2

        p = []
        keys = []

        for key in self.fp.keys():
            if not self.fp.fixed[key]:
                # we make internal parameters in case of bounding
                p.append(self.internal(self.fp[key],
                                       self.fp.lob[key], self.fp.upb[key]))
                keys.append(key)    

        res, self.fopt, self.funcalls, self.warnmess = self.minimize(chi2,
                                                        p, args=(keys,))

        self.fopt = np.sqrt(self.fopt)
        #res is internal parameters
        #restore the external parameter
        for i in range(len(p)):
            key = keys[i]
            pe = self.external(res[i], self.fp.lob[key], self.fp.upb[key])
            self.fp[key] =( pe, self.fp.lob[key],
                                self.fp.upb[key], self.fp.fixed[key])

        # update the paramscript content
        self.paramscript.script = str(self.fp)
        self._running = False
        self.fd.can_plot=True
        return

    def _callforresults(self, params, keys, funcall, iteration, fopt):
        """call back function called by minimize
        """
        for i in range(len(params)):
            key = keys[i]
            # restore external parameters in case of bounding
            pe = self.external(params[i], self.fp.lob[key], self.fp.upb[key])
            self.fp[key] = pe, self.fp.lob[key], self.fp.upb[key]

        # update parameter script
        self.funcalls = funcall
        self.iterations = iteration
        self.fopt = np.sqrt(fopt)

        return self.fittingthread.stopped   

    #===========================================================================
    # minimize
    #===========================================================================
    def minimize(self, func, x0, args=()):
        """
        """
        
        _epsilon = np.sqrt(np.finfo(float).eps)
        xtol = ftol = self.tolerance
        maxiter = self.maxiter
        maxfun = self.maxfuncalls
        callback = self._callforresults
        
        #allvecs = None
        
        if self.method == 'SIMPLEX':
            res = optimize.fmin(func, x0, args=args, xtol=xtol,
                                ftol=ftol, maxiter=maxiter,
                                maxfun=maxfun, full_output=True, disp=True,
                                retall=True, callback=callback)
            return res[0], res[1], res[3], res[4]

#        elif self.method == 'BFGS':
#            res = optimize.fmin_bfgs(func, x0, args=args, gtol=ftol, 
# norm=np.Inf,
#                      epsilon=_epsilon, maxiter=maxiter, 
# full_output=True, disp=True,
#                       retall=True, callback=callback)
#            #res = (xopt, {fopt, gopt, Hopt, func_calls, grad_calls, 
# warnflag}, <allvecs>)
#            return res[0]
#        
#        elif self.method == 'CG':
#            res = optimize.fmin_cg(func, x0, args=args, gtol=ftol, norm=np.Inf,
#                      epsilon=_epsilon, maxiter=maxiter, full_output=1, disp=1,
#                      retall=1, callback=callback)
#            return res[0]
        
        else:
            raise NotImplementedError, "Method not implemented"     
        

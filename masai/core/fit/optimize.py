#===============================================================================
# imports
#===============================================================================
import numpy as np
import scipy.optimize
from masai.core.fit.parameters import FitParameters

__all__ = ["optimize",]

#===========================================================================
# _fitting
#===========================================================================

#  Internal/external transformation
#  These transformations are used in the MINUIT package,
#  and described in detail
#  in the section 1.3.1 of the MINUIT User's Guide.

def optimize(func, fp0, args=(), constraints={}, method="SIMPLEX",
             maxfun=None, maxiter=1000,
             ftol=1e-8, xtol=1e-8,
             callback=None):
    """
    Optimization procedure

    :param func:
    :param x0:
    :param args:
    :param method:
    :param maxfun:
    :param maxiter:
    :param ftol:
    :param xtol:
    :param callback:
    :return:

    """
    global keys

    def internal_func(p, dat, fp, keys, *args):
        # restore external parameters
        for key in fp.keys():
            keysp = key.split('_')
            if keysp[0] in fp.expvars:
                if not "_".join(key[:-1]) in keys:
                    break
                ps = []
                for i in xrange(fp.expnumber):
                    k = keys.index("%s_exp%d"%(key, i))
                    ps.append(p[k])
                fp.to_external(key, ps)
            else:
                if not key in keys:
                    break
                k = keys.index(key)
                fp.to_external(key, p[k])
        return func(fp, dat, *args)

    def internal_callback(*args):
        if callback is None:
            return
        global keys
        p = args[0]
        fp = fp0.copy()
        for i in range(len(p)):
            key = keys[i]
            # restore external parameters in case of bounding
            fp.to_external(key, p[i])
            #print(fp[key], p[i])
        return callback(fp)

    if not isinstance(fp0, FitParameters):
        raise TypeError('fp0 is not of FitParameter type')

    # make internal parameters
    par = []
    keys= []
    for key in sorted(fp0.keys()):
        if not fp0.fixed[key]:
            # we make internal parameters in case of bounding
            # We also take care of the multiple experiments
            keysp = key.split('_')[0]
            if keysp in fp0.expvars:
                for i in xrange(fp0.expnumber):
                    par.append(fp0.to_internal(key, i))
                    keys.append("%s_exp%d"%(key,i))
            else:
                par.append(fp0.to_internal(key))
                keys.append(key)

    args = list(args)
    args.append(fp0)
    args.append(keys)
    if constraints:
        args.append(constraints)

    if not maxfun:
        maxfun = 4*maxiter
    if method.upper() == "SIMPLEX":
        result = scipy.optimize.fmin(internal_func, par,
                                     args=tuple(args),
                                     maxfun=maxfun, maxiter=maxiter,
                                     ftol=ftol, xtol=xtol,
                                     full_output=True, disp=False,
                                     callback=internal_callback)
        res, fopt, iterations, funcalls, warnmess = result
        # print(iterations, funcalls, warnmess)

    elif method.upper() == "HOPPING":
        result = scipy.optimize.basinhopping(internal_func, par,
                                             niter=100, T=1.0, stepsize=0.5,
                                             minimizer_kwargs={'args': tuple(args)},
                                             take_step=None,
                                             accept_test=None, callback=internal_callback,
                                             interval=50, disp=False, niter_success=None)

        # fmin(func, par, args=args, maxfun=maxfun, maxiter=maxiter, ftol=ftol, xtol=xtol,
        #                                                full_output=True, disp=False, callback=callback)
        res, fopt, warnmess = result.x, result.fun, result.message


    elif method == "XXXX":
        raise NotImplementedError("method: %s" % method)
        #TODO: implement other algorithms
    else:
        raise NotImplementedError("method: %s"%method)

    #restore the external parameter
    for i, key in enumerate(keys):
        fp0.to_external(key, res[i])

    return fp0, fopt

if __name__=="__main__":

    pass


#===============================================================================
# imports
#===============================================================================
import numpy as np
import scipy.optimize

__all__ = ["optimize", "external", "internal"]

#===========================================================================
# _fitting
#===========================================================================

#  Internal/external transformation
#  These transformations are used in the MINUIT package,
#  and described in detail
#  in the section 1.3.1 of the MINUIT User's Guide.

def optimize(func, x0, args=(), method="SIMPLEX", maxfun=5000, maxiter=5000,
             ftol=1e-8, xtol=1e-8, callback=None):
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
    par = []
    for x in x0:
        if not isinstance(x, list):
            par.append(x)
        else:
            # we make internal parameters in case of bounding
            par.append(internal(*x))

    if method.upper() == "SIMPLEX":
        result = scipy.optimize.fmin(func, par,
                                     args=args, maxfun=maxfun, maxiter=maxiter,
                                     ftol=ftol, xtol=xtol,
                                     full_output=True, disp=False,
                                     callback=callback)
        res, fopt, iterations, funcalls, warnmess = result

    elif method.upper() == "HOPPING":
        result = scipy.optimize.basinhopping(func, par,
                                             niter=100, T=1.0, stepsize=0.5,
                                             minimizer_kwargs={'args': args}, take_step=None,
                                             accept_test=None, callback=callback,
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
    xres = []
    for i, x in enumerate(x0):
        if not isinstance(x, list):
            xres.append(res[i])
        else:
            # we make internal parameters in case of bounding
            xres.append(external(res[i], x[1], x[2]))

    return xres, fopt

def internal(pe, lob=None, upb=None):
    is_lob = lob is not None # (lob > -0.1/sys.float_info.epsilon)
    is_upb = upb is not None #(upb < +0.1/sys.float_info.epsilon)

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

def external(pi, lob=None, upb=None):
    is_lob = lob is not None #(lob > -0.1/sys.float_info.epsilon)
    is_upb = upb is not None #(upb < +0.1/sys.float_info.epsilon)

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
    return pe

if __name__=="__main__":

    p = [17.184196543405619, -180.0, 180.0]
    print internal(*p)


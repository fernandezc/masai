# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.plots.plot2d_plugin
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
Define plot2d plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  Python imports
# ===============================================================================
import numpy as np
# ===============================================================================
#  Enthought imports
# ===============================================================================
from traits.api import (HasTraits, Bool, Str, Any, provides)

# ===============================================================================
# matplotlib imports
# ===============================================================================
from matplotlib import pyplot as pl
# from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ===============================================================================
# other imports
# ===============================================================================
from masai.plugins import IPlugin
import traits.has_traits
traits.has_traits.CHECK_INTERFACES = 2

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Plot2d'


# ===================================================================================
# Plot1d plugin
# ===================================================================================
@provides(IPlugin)
class Plot2d(HasTraits):

    # minimal plugin attribute -----------------------------------------------------
    description = Str("Plot2d (type: plots), output a matplotlib 2D contour plot.")
    type = Str("plot")

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
    def help(self):
        return self.__call__.__doc__

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """Plot a 2D source

        Parameters
        ----------
        fd: 2D source (such as NMR or IR source)
        ax: Axe instance, optional
            the axe where to plot
        figsize: tuple, optional, default is (6.8,6.8)
            figure size
        xlim: tuple, optional
            limit on the horizontal axis
        zlim or ylim: tuple, optional
            limit on the vertical axis
        alpha: float, optional
            alpha transparency
        linewidth or lw: float, optional
            line width
        colormap or cmap: str, optional
            color map definition
        xlabel: str, optional
            label on the horizontal axis
        ylabel: str, optional
            label on the vertical axis
        showz: bool, optional, default=True
            should we show the vertical axis
        hold: bool, optional, default is False
            if True, hold it until a new plot is performed
        """

        if args:
            # look in args
            fd = args[0]
        else:
            raise ValueError("no valid source specified\nUsage: Plot2d(fd)")

        if not fd.is_2d:
            print("This is not a 2D source.\nCannot Plot!")
            return

        # data
        x = fd.data.columns.values
        y = fd.data.index.values
        spec = fd.real

        # setup figure and axes
        figsize = kargs.get('figsize', (6.8, 6.8))

        # get the proj keyword
        proj = kargs.get('proj')

        if 'ax' not in kargs.keys():
            fig = pl.figure(figsize=figsize)
            ax = fig.add_subplot(111)
            if proj:
                # create new axes on the right and on the top of the current axes
                # The first argument of the new_vertical(new_horizontal) method is
                # the height (width) of the axes to be created in inches.
                divider = make_axes_locatable(ax)
                # print divider.append_axes.__doc__
                axx = divider.append_axes("top", 1.1, pad=0.1, sharex=ax, frameon=0, yticks=[])
                axy = divider.append_axes("right", 1.1, pad=0.1, sharey=ax, frameon=0, xticks=[])
                axx.tick_params(bottom='off', top='off')
                axy.tick_params(right='off', left='off')
                pl.setp(axx.get_xticklabels() + axx.get_yticklabels()
                        + axy.get_xticklabels() + axy.get_yticklabels(),
                        visible=False)
            xlim = [x[0], x[-1]]
            ylim = (y[0], y[-1])
            cl = self.clevels(spec.values, **kargs)

        else:
            ax = kargs.get('ax')
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            cl = ax.get_c


        # limits of the spectra
        xlim = list(kargs.get('xlim', xlim))
        xlim.sort()
        if kargs.get('xreverse', False):
            xlim.reverse()
        ylim = list(kargs.get('ylim', ylim))
        ylim.sort()
        if kargs.get('yreverse', False):
            xlim.reverse()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # attribute
        colormap = kargs.get('colormap', kargs.get('cmap','jet'))
        lw = kargs.get('linewidth', kargs.get('lw', 1.))
        alpha = kargs.get('alpha', kargs.get('alpha', .8))

        c = ax.contour(x, y, spec.values, cl, linewidths=lw, alpha=alpha)
        c.set_cmap(colormap)

        # label
        xlabel = kargs.get('xlabel', ax.get_xlabel())
        if not xlabel:
            xlabel = r'$\mathrm{time\, (%s)}$' % fd.par.units \
                if fd.datatype in ['FID', 'SER'] \
                else r'$\mathrm{\delta\, (%s)}$' % fd.par.units
        ax.set_xlabel(xlabel)

        zlabel = kargs.get('zlabel', ax.get_ylabel())
        if not zlabel:
            zlabel = r'$\mathrm{intensity\,(a.u.)}$'
        if kargs.get('showz', True):
            ax.set_ylabel(zlabel)
        else:
            ax.set_yticks([])

        # Layout
        pl.tight_layout()

        # axcommands
        commands = kargs.get('axcommands', [])
        if commands:
            for command in commands:
                com, val = command.split('(')
                val = val.split(')')[0].split(',')
                args = []
                kargs = {}
                for item in val:
                    if '=' in item:
                        k, v = item.split('=')
                        kargs[k.strip()] = eval(v)
                    else:
                        args.append(eval(item))
                print(args)
                print(kargs)
                getattr(ax, com)(*args, **kargs)

        # should be after all plot commands
        savename = kargs.get('savefig', False)
        if savename:
            pl.savefig(savename)

        if not kargs.get('hold', False):
            pl.show()

        return ax

   #===========================================================================
    # clevels
    #===========================================================================
    def clevels(self, data, **kargs):
        """Utility function to determine contours levels
        """
        # contours
        maximum = data.max()
        minimum = -maximum
        spacing = kargs.get('spacing',20)
        exponent = kargs.get('exponent',1.2)
        start = kargs.get('start',-1)

        if (exponent-1.00)<.005:
            clevelc = np.linspace(minimum, maximum, spacing)
            clevelc[clevelc.size/2-1:clevelc.size/2+1]=np.NaN
            return clevelc

        maximum = data.max()
        ms = maximum/spacing
        for xi in range(100):
            if ms*exponent**xi>maximum:
                xl = xi
                break
        if start!=0:
            clevelc = [float(start)*ms*exponent**xi
                      for xi in range(xl)]+[ms*exponent**xi for xi in range(xl)]
        else:
            clevelc = [ms*exponent**xi for xi in range(xl)]
        return clevelc

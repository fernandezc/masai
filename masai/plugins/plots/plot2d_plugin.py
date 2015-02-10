# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.plots.plot2d_plugin
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
Define plot2d plugin functionalities
"""
from __future__ import print_function, division

# ===============================================================================
#  Python imports
# ===============================================================================
import numpy as np
import sys

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
from matplotlib.ticker import MaxNLocator

# ===============================================================================
# other imports
# ===============================================================================
from masai.plugins import IPlugin
import traits.has_traits

traits.has_traits.CHECK_INTERFACES = 2
from masai.core.utils import htmldoc

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

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """
        Plot2d
        Plot a 2D source

        Parameters
        ----------
        source: 2D source (such as NMR or IR source)
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
        nc = int, default=20
            number of contour levels (between -max and +max)
            where max is the maximum amplitude of the spectra.
            if there is no negative values,
            then at the maximum nc/2 contours will be displayed
        exponent = float, default=1.2
            exponent
        start = Integer, default=-1
            start
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
            self.source = source = args[0]
        else:
            return

        if not source.is_2d:
            print("This is not a 2D source.\nCannot Plot! Use Plot1d instead", file=sys.stderr)
            return

        # data
        x = source.data.columns.values
        y = source.data.index.values
        spec = source.real

        # setup figure and axes
        figsize = kargs.get('figsize', (3.4, 3.4))
        dpi = kargs.get('dpi', 150)
        fontsize = kargs.get('fontsize', 10)

        # when using matplotlib inline
        # dpi is the savefig.dpi so we should set it here
        pl.rcParams['figure.figsize'] = figsize
        pl.rcParams['figure.dpi'] = dpi
        pl.rcParams['savefig.dpi'] = dpi
        pl.rcParams['font.size'] = fontsize
        pl.rcParams['legend.fontsize'] = int(fontsize * .8)
        pl.rcParams['xtick.labelsize'] = fontsize
        pl.rcParams['ytick.labelsize'] = fontsize

        # get the proj keyword
        proj = kargs.get('proj')
        xproj = kargs.get('xproj')
        yproj = kargs.get('yproj')

        if 'axe' not in kargs.keys():
            fig = pl.figure(figsize=figsize, dpi=dpi)
            self.fig = fig
            axe = fig.add_subplot(1, 1, 1)
            if proj or xproj or yproj:
                # create new axes on the right and on the top of the current axes
                # The first argument of the new_vertical(new_horizontal) method is
                # the height (width) of the axes to be created in inches.
                divider = make_axes_locatable(axe)
                # print divider.append_axes.__doc__
                if proj or xproj:
                    axex = divider.append_axes("top", 1.01, pad=0.01, sharex=axe, frameon=0, yticks=[])
                    axex.tick_params(bottom='off', top='off')
                    pl.setp(axex.get_xticklabels() + axex.get_yticklabels(), visible=False)
                    self.axex = axex
                if proj or yproj:
                    axey = divider.append_axes("right", 1.01, pad=0.01, sharey=axe, frameon=0, xticks=[])
                    axey.tick_params(right='off', left='off')
                    pl.setp(axey.get_xticklabels() + axey.get_yticklabels(), visible=False)
                    self.axey = axey

            xlim = [x[0], x[-1]]
            ylim = (y[0], y[-1])
            cl = self.clevels(spec.values, **kargs)
            self.axe = axe

        else:
            axe = kargs.get('axe')
            xlim = axe.get_xlim()
            ylim = axe.get_ylim()
            cl = axe.get_c

        # limits of the spectra
        xlim = list(kargs.get('xlim', xlim))
        xlim.sort()
        if hasattr(source, 'par'):
            reverse = source.par.isfreq
        else:
            reverse = False
        if kargs.get('xreverse', reverse):
            xlim.reverse()
        ylim = list(kargs.get('ylim', ylim))
        ylim.sort()
        if hasattr(source, 'par2'):
            reverse = source.par2.isfreq
        else:
            reverse = False
        if kargs.get('yreverse', reverse):
            ylim.reverse()
        axe.set_xlim(xlim)
        axe.set_ylim(ylim)

        # attribute
        colormap = kargs.get('colormap', kargs.get('cmap', 'jet'))
        lw = kargs.get('linewidth', kargs.get('lw', 1.))
        alpha = kargs.get('alpha', kargs.get('alpha', .8))

        c = axe.contour(x, y, spec.values, cl, linewidths=lw, alpha=alpha)
        c.set_cmap(colormap)

        # label
        xlabel = kargs.get('xlabel', axe.get_xlabel())
        if hasattr(source, 'par'):
            ppm = source.par.isfreq
        if not xlabel:
            xlabel = r'$\mathrm{time\, (%s)}$' % source.par.units \
                if source.datatype in ['FID', 'SER'] and not ppm \
                else r'$\mathrm{\delta\, (%s)}$' % source.par.units
        axe.set_xlabel(xlabel)

        #
        zlabel = kargs.get('zlabel', axe.get_ylabel())
        if hasattr(source, 'par2'):
            ppm = source.par2.isfreq
        if not zlabel and hasattr(source, 'par2'):
            zlabel = r'$\mathrm{time\, (%s)}$' % source.par2.units \
                if source.datatype in ['SER'] and not ppm \
                else r'$\mathrm{\delta\, (%s)}$' % source.par2.units
        if kargs.get('showz', True):
            axe.set_ylabel(zlabel)
        else:
            axe.set_yticks([])

        axe.xaxis.set_major_locator(MaxNLocator(5))
        axe.yaxis.set_major_locator(MaxNLocator(5))

        if proj == 'sky' or xproj == 'sky':
            dxp = spec.values.max(axis=0)
            xproj, = axex.plot(x, dxp)
        if proj == 'sky' or yproj == 'sky':
            dyp = spec.values.max(axis=1)
            yproj, = axey.plot(dyp, y)
        if proj == 'sum' or 'sum' == xproj:
            dxp = spec.values.sum(axis=0)
            dxp = dxp / (dxp.max() - dxp.min())
            xproj, = axex.plot(x, dxp)
        if 'sum' == proj or 'sum' == yproj:
            dyp = spec.values.sum(axis=1)
            dyp = dyp / (dyp.max() - dyp.min())
            yproj, = axey.plot(dyp, y)

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
                getattr(axe, com)(*args, **kargs)

        if not kargs.get('hold', False):
            self.show()

        # should be after all plot commands
        savename = kargs.get('savename', None)
        if savename:
            self.savename = savename
            self.fig.savefig(savename)

    def show(self):
        pl.show()

    def savefig(self, savename=None):
        if savename is None:
            savename = self.savename
        self.fig.savefig(savename)
        self.savename = savename

        #===========================================================================

    # clevels
    #===========================================================================
    def clevels(self, data, **kargs):
        """Utility function to determine contours levels
        """
        # contours
        maximum = data.max()
        minimum = -maximum
        if kargs.has_key('spacing'):
            print("'spacing' is deprecated, use 'nc' instead.", file=sys.stderr)
        spacing = kargs.get('nc',kargs.get('spacing', 20))
        exponent = kargs.get('exponent', 1.2)
        start = kargs.get('start', -1)

        if (exponent - 1.00) < .005:
            clevelc = np.linspace(minimum, maximum, spacing)
            clevelc[clevelc.size / 2 - 1:clevelc.size / 2 + 1] = np.NaN
            return clevelc

        maximum = data.max()
        ms = maximum / spacing
        for xi in range(100):
            if ms * exponent ** xi > maximum:
                xl = xi
                break
        if start != 0:
            clevelc = [float(start) * ms * exponent ** xi
                       for xi in range(xl)] + [ms * exponent ** xi for xi in range(xl)]
        else:
            clevelc = [ms * exponent ** xi for xi in range(xl)]
        return clevelc

    def _repr_html_(self):
        if not self.axe:
            return htmldoc(self.__init__.__doc__)
        else:
            return

if __name__ == '__main__':
    """
    test plugin
    """
    pass



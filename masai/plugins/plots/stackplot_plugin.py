# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.plots.stackplot_plugin
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
Define plot1d plugin functionalities
"""
from __future__ import print_function, division

#
# Python imports
#
from IPython.display import display, display_pretty
import numpy as np
import pylab as pl
from traits.api import (HasTraits, Bool, Str, Int, Any, provides)
from pandas import DataFrame, Series, concat
import matplotlib as mpl
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.ticker import MaxNLocator
from masai.core.source import Source

# ===============================================================================
# other imports
# ===============================================================================
from masai.plugins import IPlugin
import traits.has_traits
from masai.core.utils import htmldoc
traits.has_traits.CHECK_INTERFACES = 2

DEBUG = False

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Stackplot'

all = ['Stackplot', ]

# ===================================================================================
# Stackplot plugin
# ===================================================================================
@provides(IPlugin)
class Stackplot(HasTraits):
    # minimal plugin attribute -----------------------------------------------------
    description = Str("Stackplot (type: plots), output several stacked matplotlib 1D plot.")
    type = Str("plot")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)
    savename = Str
    dpi = Int

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    def _repr_html_(self):
        if not self.source:
            return htmldoc(self.__init__.__doc__)
        else:
            return self.message

    def __repr__(self):
        return self.__init__.__doc__

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """ Stackplot

        Stack plot of 1D sources

        Parameters
        ----------
        sources: list of 1D sources (such as NMR or IR source), DataFrame, Series or a 2D source

        figure: Figure instance, optioanl
        axe: Axe instance, optional
            the axe where to plot
        figsize: tuple, optional, default is (6.8,3.4)
            figure size
        dpi: int, optional, default=150
            dot per inch
        fontsize: int, optional, default=10
            font size
        imag: bool, optional, default False
            By default real part is shown. Set to True to display the imaginary part
        xlim: tuple, optional
            limit on the horizontal axis
        zlim or ylim: tuple, optional
            limit on the vertical axis
        color or c: matplotlib valid color or str, optional
            color of the line or 'map'
        colormap: colormap
            if color==map, choice of the color map to use
        normalize: bool, optional, default None
            normalize
        colorbar: boolean, optional, default=True
            whether or not to display the colorbar
        alpha: float, optional, default 0.5
            opacity between 0 and 1
        linewidth or lw: float, optional
            line width
        linestyle or ls: str, optional
            line style definition
        xlabel: str, optional
            label on the horizontal axis
        zlabel or ylabel: str, optional
            label on the vertical axis
        showz: bool, optional, default=True
            should we show the vertical axis
                colorbar = kargs.get('colorbar', True)
        fill: bool, optional, default=True
            fill below curves
        offset: float, optional, default=0
            separation between curves
        hold: bool, optional, default is False
            if True, hold it until a new plot is performed
        """
        if args:
            # look in args
            self.sources = sources = args[0]
        else:
            return

        # in principle, if we want a stack plot we must provide a list of source or a 2D spectrum
        if isinstance(sources, list):
            # ok we have a list
            # is-it a list of source object?
            if isinstance(sources[0], Source):
                listdata = [source.data for source in sources]
            # or of pandas DataFrame?
            elif isinstance(sources[0],DataFrame):
                listdata = sources
            # or of pandas Series?
            elif isinstance(sources[0], Series):
                # this routine needs DataFrame input, not Series
                listdata = [DataFrame(source).T for source in sources]
            else:
                # cannot handle object which are not pandas object.
                # Todo: may be we could provide the x scale to allow numpy arrays in general
                raise TypeError('The input must be a list of pandas DataFrame or Series')
            # ok now we have list of dataframes
            # lets create a 2D dataframe
            data = concat(listdata)
            # todo: set the index
            y = kargs.get('index')

        elif isinstance(sources, Source):
            # it should be a 2D source
            if sources.is_2d:
                data = sources.data
                y = data.index.values
            else:
                raise TypeError("I can't interpret the type of the input")

        if y is None:
            y = np.arange(data.shape[0])
        # setup figure and axes
        figsize = kargs.get('figsize', (3.4, 3.4))
        dpi = kargs.get('dpi', 150)
        fontsize = kargs.get('fontsize', 10)

        # when using matplotlib inline
        # dpi is the savefig.dpi so we should set it here
        pl.rcParams['figure.figsize']=figsize
        pl.rcParams['figure.dpi']=dpi
        pl.rcParams['savefig.dpi']=dpi
        pl.rcParams['font.size']=fontsize
        pl.rcParams['legend.fontsize'] = int(fontsize*.8)
        pl.rcParams['xtick.labelsize']=fontsize
        pl.rcParams['ytick.labelsize']=fontsize

        if 'axe' not in kargs.keys():
            fig = pl.figure(figsize=figsize, dpi=dpi)
            axe = fig.add_subplot(1, 1, 1)
        else:
            axe = kargs.get('axe')

        # attribute
        c = kargs.get('color', kargs.get('c', 'map'))
        alpha = kargs.get('alpha', 0.5)
        lw = kargs.get('linewidth',  kargs.get('lw', 1.))
        ls = kargs.get('linestyle', kargs.get('lw', 'solid'))
        colormap = kargs.get('colormap', kargs.get('cm', 'jet'))
        normalize = kargs.get('normalize', None)
        colorbar = kargs.get('colorbar', True)
        fill = kargs.get('fill', True)
        offset = kargs.get('offset', 0.0)

        # labels
        xlabel = kargs.get('xlabel', axe.get_xlabel())
        if not xlabel and hasattr(sources, 'par'):
            xlabel = sources.par.label
        elif not xlabel and hasattr(sources[0], 'par'):
            xlabel = sources[0].par.label
        axe.set_xlabel(xlabel)

        ylabel = kargs.get('ylabel', None)
        if not ylabel and hasattr(sources, 'par2'):
            ylabel = sources.par2.label

        zlabel = kargs.get('zlabel', axe.get_ylabel())
        if not zlabel:
            zlabel = r'$\mathrm{intensity\,(a.u.)}$'

        spec = data.values.real
        x = data.columns.values

        yoffset = spec.max()*offset/100.
        speco = (spec.T + np.arange(spec.shape[0])*yoffset).T

        if c != 'map':
            # plot the source data
            lines = axe.plot(x, speco.T)
            for line in lines:
               line.set_color(c)
               line.set_linewidth(lw)
               line.set_linestyle(ls)

        else: # map
            norm = mpl.colors.Normalize(vmin=y[0], vmax=y[-1])  # we normalize to the max time

            #TODO: add a kargs
            if normalize:
                norm.vmax = normalize[-1]
                norm.vmin = normalize[0]

            if fill:
                fun = PolyCollection
                lw /= 5.
            else:
                fun = LineCollection

            line_segments = fun([list(zip(x, speco[i].tolist())) for i in xrange(len(y))][::-1],
                                        linewidths=lw,
                                        linestyles=ls,
                                        alpha=alpha,
                             )

            line_segments.set_array(y[::-1])
            line_segments.set_cmap(colormap)
            line_segments.set_norm(norm)

            axe.add_collection(line_segments)

            if colorbar:
                fig = pl.gcf()
                axcb = fig.colorbar(line_segments, ax=axe)
                if ylabel:
                    axcb.set_label(ylabel)


        # limits of the spectra
        if not axe in kargs:
            xlim = (x[0], x[-1])
            mi,ma = spec.min(),spec.max()
            amp = ma-mi
            zlim = (mi-amp*.05,ma+amp*.05)
        else:
            xlim = axe.get_xlim()
            zlim = axe.get_ylim()

        xlim = list(kargs.get('xlim', xlim))
        xlim.sort()
        xl = [x[0], x[-1]]
        xl.sort()
        xlim[-1] = min(xlim[-1], xl[-1])
        xlim[0] = max(xlim[0], xl[0])
        if hasattr(sources, 'par'):
            reverse = sources.par.isfreq
        else:
            reverse = False
        if kargs.get('reverse', reverse):
            xlim.reverse()
        zlim = list(kargs.get('zlim', kargs.get('ylim', zlim)))
        zlim.sort()
        axe.set_xlim(xlim)
        axe.set_ylim(zlim)
        axe.xaxis.set_major_locator(MaxNLocator(4))
        axe.yaxis.set_major_locator(MaxNLocator(4))

        if kargs.get('showz', True):
            axe.set_ylabel(zlabel)
        else:
            axe.set_yticks([])

        # Layout
        pl.tight_layout()

        # commands
        commands = kargs.get('commands', [])
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
                if DEBUG:
                    print(args)
                    print(kargs)
                getattr(axe, com)(*args, **kargs)

        # should be after all plot commands
        savename = kargs.get('savefig', False)
        if savename:
            self.savefig(savename)

        if not kargs.get('hold', False):
            self.show()

        self.axe = axe

    def show(self):
        pl.show()

    def savefig(self, savename=None):
        if savename is None:
            savename = self.savename
        pl.savefig(savename)
        self.savename = savename


if __name__ == '__main__':
    """
    test plot1d plugin

    """
    # need this to test  directly
    import os
    from masai.plugins.list import Bruker, Process

    Bruker = Bruker()  # necessary for direct call
    Stackplot = Stackplot()
    Process = Process()


    #DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    #user = 'christian'
    #name = 'example1'
    #expno = 10
    #path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    #fd = Bruker(path)

    #sources = [fd, fd, fd, fd]
    ## plot a Source
    #ax = Stackplot(sources, poly=False)

    # plot a Dataframe
    #ax = Stackplot(sources, color='r')

    DATADIR = os.path.expanduser(os.path.join('~','Dropbox','storage','bruker','data'))


    sources = Bruker(
                data_dir=DATADIR,
                user = 'siblanixenon',
                name = 'NaY1-QIN',
                expnos = [25,49,65,85,105,124,144,164,184,194,209,224,242],
                varpars = "TE",
                )

    Process(sources, transform='em 50; zf 0; ft; cl -100 300; pk --auto')# --fit_phc1 --byrow')

    Stackplot(sources, offset=5., ylabel='temperature (K)', colorbar=False)



# -*- coding: utf-8 -*-
#
# ===============================================================================
# masai.plugins.plots.plot1d_plugin
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
import matplotlib.pyplot as pl

from traits.api import (HasTraits, Bool, Str, Int, Any, provides)
from pandas import DataFrame, Series
from matplotlib.ticker import MaxNLocator

#
# other imports
#
from masai.plugins import IPlugin
from masai.core.utils import htmldoc, is_sequence

import traits.has_traits

traits.has_traits.CHECK_INTERFACES = 2

DEBUG = False

# ===============================================================================
# Plugin signature
# ===============================================================================
__plugin_main__ = 'Plot1d'

# ===================================================================================
# Plot1d plugin
# ===================================================================================
@provides(IPlugin)
class Plot1d(HasTraits):
    # minimal plugin attribute -----------------------------------------------------
    plugin_info = Str("Plot1d (type: plots), output a matplotlib 1D plot.")
    plugin_type = Str("plot")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)
    savename = Str
    dpi = Int

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    def _repr_html_(self):
        if not self.axe:
            return htmldoc(self.__init__.__doc__)
        else:
            return

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __init__(self, *args, **kargs):
        """
        Plot1d
        Plot a 1D source
        Usage: Plot1d(source, **kargs)

        Parameters
        ----------
        source: source(s)
            A source or a list of source (such as NMR or IR source)
        axe: Axe instance, optional
            the axe where to plot
        subplot: str
            a subplot position e.g., 121
        dpi: int, optional
            the number of pixel per inches
        figsize: tuple, optional, default is (3.4, 1.7)
            figure size
        fontsize: int, optional
            font size in pixels, default is 10
        imag: bool, optional, default False
            By default real part is shown. Set to True to display the imaginary part
        xlim: tuple, optional
            limit on the horizontal axis
        zlim or ylim: tuple, optional
            limit on the vertical axis
        color or c: matplotlib valid color, optional
            color of the line
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
        hold: bool, optional, default is False
            if True, hold it until a new plot is performed
        plot_model:Bool,
            plot model data if available
        modellinestyle or modls: str,
            line style of the model
        offset: float,
            offset of the model individual lines
        commands: str,
            matplotlib commands to be executed
        show_zero: boolean, optional
            show the zero basis
        savename: str,
            name of the file to save the figure
        savefig: Bool,
            save the fig if savename is defined, should be executed after all other commands
        vshift: float, optional
            vertically shift the line from its baseline
        """

        multiplots = False

        if args:
            # look in args
            source = args[0]
            if isinstance(source, DataFrame):
                multiplots = False
            elif is_sequence(source):
                # we will draw subplots
                multiplots = True
        else:
            # we will display the help
            return

        imag = kargs.get('imag', False)
        specs = []

        if not multiplots:

            if isinstance(source, DataFrame):
                data = source

            elif isinstance(source, Series):
                # this routine needs DataFrame input, not Series
                data = DataFrame(source).T

            else:
                data = source.data

            spec = data.values.real if not imag else data.values.imag

            if (spec.ndim >= 2) & (spec.shape[0] > 1):
                # not well suitable for this
                print('Cannot plot a 2D source.\n Use Plot2d.')
                return
            specs = [spec]

        else:

            for s in source:
                if isinstance(s, DataFrame):
                    data = s
                elif isinstance(s, Series):
                    # this routine needs DataFrame input, not Series
                    data = DataFrame(s).T
                else:
                    data = s.data
                spec = data.values.real if not imag else data.values.imag
                specs.append(spec)

        # setup figure and axes
        figsize = kargs.get('figsize', (3.4, 1.7))
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

        if 'axe' not in kargs.keys():
            fig = self.fig = pl.figure(figsize=figsize, dpi=dpi)
            axe = fig.add_subplot(1, 1, 1)
        else:
            axe = kargs.get('axe')

        # data
        x = data.columns.values

        vshift = kargs.get('vshift', 0.0)

        # plot the source data
        for i, spec in enumerate(specs):
            line, = axe.plot(x, spec.T + vshift*float(i+1))

        # attribute
        c = kargs.get('color', kargs.get('c'))
        if c:
            line.set_color(c)
        lw = kargs.get('linewidth', kargs.get('lw', 1.))
        if lw:
            line.set_linewidth(lw)
        ls = kargs.get('linestyle', kargs.get('ls', '-'))
        if ls:
            line.set_linestyle(ls)

        # plot the models
        offset = 0.

        if not multiplots and kargs.get('plot_model', False) and source.modeldata is not None:
                # we now plot the model if it exists
                mdata = source.modeldata
                offset = (mdata.values.max()) * (float(kargs.get('offset', 1.0)) / 100.)
                # print('offset',offset)
                mls = kargs.get('modellinestyle', kargs.get('modls', '-'))
                for model in mdata.index:
                    l, = axe.plot(x, mdata.loc[model], alpha=.75)
                    l.set_linestyle(mls)
                    l.set_linewidth(lw * .8)
                    if model == 'baseline':
                        l.set_color('k')
                        l.set_linestyle(':')
                        l.set_linewidth(lw * .4)
                        l.set_ydata(l.get_ydata() - offset)
                    elif model != 'modelsum':
                        l.set_linestyle('-')
                        l.set_linewidth(lw * .4)
                        l.set_ydata(l.get_ydata() - offset)


        # limits of the spectra
        xlim = list(kargs.get('xlim', axe.get_xlim()))
        xlim.sort()
        xl = [x[0], x[-1]]
        xl.sort()
        xlim[-1] = min(xlim[-1], xl[-1])
        xlim[0] = max(xlim[0], xl[0])
        if hasattr(source, 'par'):
            reverse = source.par.isfreq
        else:
            reverse = False
        if kargs.get('reverse', reverse):
            xlim.reverse()

        zlim = list(kargs.get('zlim', kargs.get('ylim', axe.get_ylim())))
        zlim.sort()
        if offset > 0:
            zlim[0] = zlim[0] - offset
        axe.set_xlim(xlim)
        axe.set_ylim(zlim)
        axe.xaxis.set_major_locator(MaxNLocator(5))
        axe.yaxis.set_major_locator(MaxNLocator(5))

        # label
        xlabel = kargs.get('xlabel', axe.get_xlabel())
        if not xlabel and hasattr(source, 'par'):
            xlabel = source.par.label
        axe.set_xlabel(xlabel)

        zlabel = kargs.get('zlabel', kargs.get('ylabel', axe.get_ylabel()))
        if not zlabel:
            zlabel = r'$\mathrm{intensity\,(a.u.)}$'

        if kargs.get('showz', True):
            axe.set_ylabel(zlabel)
        else:
            axe.set_yticks([])

        if kargs.get('show_zero', False):
            axe.haxlines()

        # Layout
        #pl.tight_layout()

        # commands
        commands = kargs.get('commands', [])
        if commands:
            for command in commands:
                com, val = command.split('(')
                val = val.split(')')[0].split(',')
                ags = []
                kws = {}
                for item in val:
                    if '=' in item:
                        k, v = item.split('=')
                        kws[k.strip()] = eval(v)
                    else:
                        ags.append(eval(item))
                if DEBUG:
                    print(ags)
                    print(kws)
                getattr(axe, com)(*ags, **kws)

        # should be after all plot commands
        savename = kargs.get('savefig', False)
        if savename:
            self.savefig(savename)

        if kargs.get('display', True):
            if not kargs.get('hold', False):
                    self.show()

        self.axe = axe

        savename = kargs.get('savename', None)
        if savename:
            self.savename = savename
            self.fig.savefig(savename)

    def show(self):
        pl.show()

    def savefig(self, savename=None):
        if savename is None:
            savename = self.savename
        try:
            self.fig.savefig(savename)
        except AttributeError:
            # may be the fig is not known
            self.fig = self.axe.get_figure()
            self.fig.savefig(savename)
        self.savename = savename


if __name__ == '__main__':

    from masai.api import *
    DATADIR = '/Users/christian/Dropbox/storage/bruker/data'
    name = 'RUB-17-2015'
    user = 'anaP'

    s=[]
    for expno in [1,2,3,60,62]:
        s.append(Bruker(data_dir=DATADIR, user=user, name=name, expno=expno))
        # print(s[-1].par.NS, s[-1].par.RG, s[-1].norm)
        #    Process(s[-1],
        #        transform='em 1; zf 3; ft; cl -105 -75; pk 337.167 1.938 --auto; ab --mode poly')
        Process(s[-1],
                transform='em 1; zf 3; ft; cl -110 -75; pk 337.167 1.938 --auto')

    p = Plot1d(s, xlim=(-75,-100), display=False)

    source = Bruker(data_dir=DATADIR, user = 'user', name = 'SPEC1D', expno = 10)
    Process(source, transform="em 50; zf 1; ft; pk --auto --fit_phc1 --bound_phc1 360 --pivot 4.5")
    p = Plot1d(source, hold=True, display=False)
    Process(source, "ab --mode poly")
    p = Plot1d(source.par.baseline, axe=p.axe, xlim=(30,-10), reverse=True, display=False)


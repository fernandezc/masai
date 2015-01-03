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
from IPython.display import display, display_pretty
import pylab as pl
from traits.api import (HasTraits, Bool, Str, Int, Any, provides)
from pandas import DataFrame, Series
from matplotlib.ticker import MaxNLocator

# ===============================================================================
# other imports
# ===============================================================================
from masai.plugins import IPlugin
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
    description = Str("Plot1d (type: plots), output a matplotlib 1D plot.")
    type = Str("plot")

    # specific attributes ----------------------------------------------------------

    silent = Bool(True)
    savename = Str
    dpi = Int

    _ = Any

    # *******************************************************************************
    # public functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def print_doc(self):
        return display(self.description)

    # -------------------------------------------------------------------------------
    def help(self):
        return self.__call__.__doc__

    def _repr_html_(self):
        import re
        import string
        p = re.compile("(.*:)?(.*)", re.MULTILINE)
        m = p.findall(self.__call__.__doc__)
        html = ""
        start = False
        for mg in m:
            if mg[0].strip() != "":
                html += "<b>%s</b>:%s<br/>" % mg
            else:
                s = string.join(mg," ")
                if s.strip() and start:
                    html += "&nbsp;"*4+"%s<br/>" % s
                if "----" in s:
                    start = True
                    html = "<h3>Parameters</h3>"
        return html

    # *******************************************************************************
    # Private functions
    # *******************************************************************************

    # -------------------------------------------------------------------------------
    def __call__(self, *args, **kargs):
        """ Plot a 1D source

        Parameters
        ----------
        fd: source (such as NMR or IR source)
        figure: Figure instance, optioanl
        axe: Axe instance, optional
            the axe where to plot
        figsize: tuple, optional, default is (6.8,3.4)
            figure size
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
        """
        if args:
            # look in args
            source = args[0]
        else:
            raise ValueError("no a valid input\nUsage: Plot1d(<source>)")

        imag = kargs.get('imag', False)
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

        # setup figure and axes
        figsize = kargs.get('figsize', (3.4, 1.7))
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

        # data
        x = data.columns.values

        # plot the source data
        line, = axe.plot(x, spec.T)

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
        if kargs.get('plot_model', False) and source.modeldata is not None:
            # we now plot the model if it exists
            mdata = source.modeldata
            offset = (mdata.values.max())*(float(kargs.get('offset', 1.0))/100.)
            # print('offset',offset)
            for model in mdata.index:
                l, = axe.plot(x, mdata.loc[model],alpha=.75)
                l.set_linestyle('-')
                l.set_linewidth(lw*.8)
                if model=='baseline':
                    l.set_color('k')
                    l.set_linestyle(':')
                    l.set_linewidth(lw*.4)
                    l.set_ydata(l.get_ydata()-offset)
                elif model != 'modelsum':
                    l.set_linestyle('-')
                    l.set_linewidth(lw*.4)
                    l.set_ydata(l.get_ydata()-offset)


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
        axe.set_xlim(xlim)
        axe.set_ylim(zlim)
        axe.xaxis.set_major_locator(MaxNLocator(4))
        axe.yaxis.set_major_locator(MaxNLocator(4))

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
            pl.show()

        return axe

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
    from masai.plugins.sources.bruker_plugin import Bruker

    Bruker = Bruker()  # necessary for direct call
    Plot1d = Plot1d()

    DATADIR = os.path.expanduser(os.path.join('~', 'Dropbox', 'storage', 'bruker', 'data'))
    user = 'christian'
    name = 'example1'
    expno = 10
    path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
    fd = Bruker(path)

    # plot a Source
    ax = Plot1d(fd, hold=True)
    ax = Plot1d(fd, axe=ax, imag=True, xlim=(0, 4), ylim=(-1, 1), c='orange')

    # plot a Dataframe
    #Plot1d(fd.data)

    #pl.show()

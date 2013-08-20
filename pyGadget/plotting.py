# plotting.py
# Jacob Hummel
"""
This module contains classes for various types of plots.
"""
import multiprocessing
from matplotlib import pyplot 

class Plot(object):
    """
    Base class for plots.
    """
    def __init__(self, snapshot, **fig_args):
        super(Plot,self).__init__()
        self.snapshot = snapshot
        self.figure = pyplot.figure(**fig_args)
        self.figure.clf()

    def save(self,path):
        self.figure.savefig(path,
                            bbox_extra_artists=(self.title,),
                            bbox_inches='tight')

class Phase(Plot):
    """
    Class for phase-space style plots.
    """
    def __init__(self, snapshot, figsize=(12,8), **fig_args):
        super(Phase,self).__init__(snapshot,figsize=figsize,**fig_args)
        self.axes = self.figure.add_subplot(111)
        self._set_phase_dict()
        z = snapshot.header.Redshift
        self.title = self.figure.suptitle('Redshift: %.2f' %(z,))

    def _set_phase_dict(self):
        self._phase_plots = {'temp':self.temp,
                             'electron_frac':self.electron_frac,
                             'h2frac':self.h2frac,
                             'HDfrac':self.HDfrac}

    def plot(self, key, **kwargs):
        self._phase_plots[key](self.axes, **kwargs)

    def _hexbin(self, ax, x, y, **kwargs):
        grid = kwargs.pop('gridsize',500)
        bins = kwargs.pop('bins','log')
        xscale = kwargs.pop('xscale','log')
        yscale = kwargs.pop('yscale','log')
        mincnt = kwargs.pop('mincnt',1)
        ax.hexbin(x, y, gridsize=grid, bins=bins, xscale=xscale,
                  yscale=yscale, mincnt=mincnt, **kwargs)
        return ax

    def temp(self, ax, **kwargs):
        dens = self.snapshot.gas.get_number_density()
        temp = self.snapshot.gas.get_temperature()
        if kwargs.pop('parallel', False):
            self.snapshot.gas.cleanup('ndensity','temp')
        ax = self._hexbin(ax, dens, temp, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(2e-3, 1e12)
        ax.set_ylim(10, 2e4)
        ax.set_xlabel('n [cm$^{-3}$]')
        ax.set_ylabel('Temperature [K]')
        return ax

    def electron_frac(self, ax, **kwargs):
        dens = self.snapshot.gas.get_number_density()
        efrac = self.snapshot.gas.get_electron_fraction()
        if kwargs.pop('parallel', False):
            self.snapshot.gas.cleanup('ndensity','electron_frac')
        ax = self._hexbin(ax, dens, efrac, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-2, 1e12)
        ax.set_ylim(1e-12, 1e-2)
        ax.set_xlabel('n [cm$^{-3}$]')
        ax.set_ylabel('f$_{e^-}$')
        return ax

    def h2frac(self, ax, **kwargs):
        dens = self.snapshot.gas.get_number_density()
        h2frac = self.snapshot.gas.get_H2_fraction()
        if kwargs.pop('parallel', False):
            self.snapshot.gas.cleanup('ndensity','h2frac')
        ax = self._hexbin(ax, dens, h2frac, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-2, 1e12)
        ax.set_ylim(1e-7,2)
        ax.set_xlabel('n [cm$^{-3}$]')
        ax.set_ylabel('f$_{H_2}$')
        return ax

    def HDfrac(self, ax, **kwargs):
        dens = self.snapshot.gas.get_number_density()
        HDfrac = self.snapshot.gas.get_HD_fraction()
        if kwargs.pop('parallel', False):
            self.snapshot.gas.cleanup('ndensity','HDfrac')
        ax = self._hexbin(ax, dens,HDfrac, **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-2, 1e12)
        ax.set_ylim(1e-11,1e-4)
        ax.set_xlabel('n [cm$^{-3}$]')
        ax.set_ylabel('f$_{HD}$')

class Radial(Plot):
    """
    Class for plotting various properties vs 'radius'.
    """
    def __init__(self, snapshot, **fig_args):
        super(Radial,self).__init__(snapshot, **fig_args)
        self.axes = self.figure.add_subplot(111)

class Quad(Phase, Radial):
    def __init__(self, snapshot, figsize=(12,8), **fig_args):
        super(Quad,self).__init__(snapshot,figsize=figsize,**fig_args)
        z = snapshot.header.Redshift
        self.title = self.figure.suptitle('Redshift: %.2f' %(z,))
        ax1 = self.figure.add_subplot(221)
        ax2 = self.figure.add_subplot(222)
        ax3 = self.figure.add_subplot(223)
        ax4 = self.figure.add_subplot(224)
        self.axes = [ax1,ax2,ax3,ax4]
        self._set_phase_dict()

    def plot(self, *keys, **kwargs):
        if len(keys) != 4:
            raise KeyError("Quad.plot requires exactly 4 keys!")
        for key in keys:
            if key not in self._phase_plots.keys():
                raise KeyError("Unrecognized plot request!: "+key)

        if kwargs.get('parallel', False):
            jobs = []
            for i,key in enumerate(keys):
                proc = multiprocessing.Process(target=self._phase_plots[key],
                                               args = (self.axes[i],),
                                               kwargs = kwargs)
                jobs.append(proc)
                proc.start()
            for process in jobs:
                process.join()
        else:
            for i,key in enumerate(keys):
                self._phase_plots[key](self.axes[i], **kwargs)

        density_labels = (1e-2,1e0,1e2,1e4,1e6,1e8,1e10,1e12)
        for ax in self.axes:
            ax.set_xticks(density_labels)
        self.figure.subplots_adjust(top=0.94, left=0.085, right=.915)

class Image(Plot):
    """
    Class for image-style plots.
    """
    def __init__(self, snapshot, **kwargs):
        super(Image,self).__init__(snapshot,**kwargs)



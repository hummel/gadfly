# plotting.py
# Jacob Hummel
"""
This module contains classes for various types of plots.
"""
import numpy
import multiprocessing
from matplotlib import pyplot 

import units
import analyze
import visualize

class Plot(object):
    """
    Base class for plots.
    """
    def __init__(self, snapshot, **fig_args):
        super(Plot,self).__init__()
        self.snapshot = snapshot
        self.figure = pyplot.figure(1,**fig_args)
        self.figure.clf()
        self.title = None

    def save(self,path):
        if self.title:
            self.figure.savefig(path, bbox_extra_artists=(self.title,),
                                bbox_inches='tight')
        else:
            self.figure.savefig(path, bbox_inches='tight')

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
        ax.axhline(2.725 * (snapshot.header.Redshift + 1),
                   linestyle='--', color='k')
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
        ax.set_xlim(2e-3, 1e12)
        ax.set_ylim(2e-11,1e-4)
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
    def __init__(self, snapshot, figsize=(12,9), **fig_args):
        self.track_sinks = fig_args.pop('track_sinks',False)
        super(Image,self).__init__(snapshot,figsize=figsize,**fig_args)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_axis_off()
        self.cbar = None

    def annotate_axes(self, scale=None):
        self.axes.text(.01, .96, 'z: %.2f' %self.snapshot.header.Redshift,
                        color='white', fontsize=18,
                        transform=self.axes.transAxes)
        if scale:
            boxsize = "".join(ch if ch.isdigit() else "" for ch in scale)
            unit = "".join(ch if not ch.isdigit() else "" for ch in scale)
            coordinates = self.snapshot.sim.units.coordinate_system
            text = '{} {} ({})'.format(boxsize,unit,coordinates)
            self.axes.text(.01,.01, text, color='white', fontsize=18,
                            transform=self.axes.transAxes)
        if self.track_sinks:
            if self.snapshot.sim.tsink:
                t_acc = self.snapshot.header.Time * units.Time_yr \
                    - self.snapshot.sim.tsink
                self.axes.text(.75, .96, 't$_{form}$: %.0f yr' %t_acc,
                                color='white', fontsize=18,
                                transform=self.axes.transAxes)
                for sink in self.snapshot.sinks:
                    self.axes.plot(sink.y, sink.x, 'ko', ms=5, mew=1)
                    self.axes.text(sink.y+15, sink.x+12, ' %.1f' %sink.mass)

    def density(self, scale, viewpoint, **kwargs):
        x,y,z = visualize.project(self.snapshot, 'ndensity',
                                  scale, viewpoint, **kwargs)
        ax = kwargs.pop('axis',self.axes)
        ax.cla()
        ax.set_axis_off()
        img = ax.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()],
                        cmap=pyplot.cm.RdGy_r,origin='lower')
        clim = kwargs.pop('clim',None)
        if clim:
            img.set_clim(clim[0],clim[1])

        cbar = kwargs.pop('colorbar', True)
        if cbar:
            if not self.cbar:
                self.cbar = self.figure.colorbar(img)
            if clim:
                self.cbar.set_clim(clim[0],clim[1])
                self.cbar.set_ticks(range(clim[0], clim[1]+1))
            self.cbar.set_label('Log Number Density [cm$^{-3}$]')
        self.annotate_axes(scale)
        self.axes.set_xlim(x.min(), x.max())
        self.axes.set_ylim(y.min(), y.max())
        pyplot.draw()

#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy

from scipy import weave
from scipy.weave import converters
from matplotlib import pyplot
from matplotlib import cm
from matplotlib.mlab import griddata
import pp

from pyGadget import units
from pyGadget import snapshot
from pyGadget import visualize
import statistics

if len(sys.argv) < 4:
    print 'Usage: python density_vis.py (simulation name) (beginning snapshot)'\
        ' (final snapshot)'
    sys.exit()

length_unit = units.Length_AU
pps = 1000 # 'pixels' per side
hsml_factor = 1.5
simulation = sys.argv[1]
path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'
boxsize = 2e3

#pyplot.ioff()
job_server = pp.Server()

start = int(sys.argv[2])
stop = int(sys.argv[3])+1
for i in xrange(start,stop):
    fname = path + '{:0>3}'.format(i) + '.hdf5'
    snap = snapshot.File(fname)
    redshift = snap.header.Redshift
    for suffix in ['-dens-xy.png','-dens-xz.png','-dens-yz.png']:
        wpath = write_dir + '{:0>4}'.format(i) + suffix
        view = suffix[-6:-4]
        x,y,z = visualize.density(snap, view, boxsize, 1., length_unit, 
                                  job_server, pps, hsml_factor)
        z = numpy.log10(z)
        zmin,zmax = (1e9,7e11)

        print 'Plotting',view+'...'
        fig = pyplot.figure(1,(12,12))
        fig.clf()
        pyplot.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()],
                      cmap=cm.jet)
        pyplot.clim(numpy.log10(zmin),numpy.log10(zmax))
        ax = pyplot.gca()
        for sink in snap.sinks:
            ax.plot(sink[1], -sink[0], 'ko') # 90-degree rotation
        ax.set_xlim(x.min(),x.max())
        ax.set_ylim(y.min(),y.max())
        ax.set_xlabel('AU')
        ax.set_ylabel('AU')
        ax.set_title('Redshift: %.7f' %redshift)
        pyplot.draw()
        pyplot.savefig(wpath, 
                       bbox_inches='tight')
    snap.close()
job_server.destroy()


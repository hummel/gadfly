#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy

from scipy import weave
from scipy.weave import converters
from matplotlib import pyplot
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
boxsize = 1e3

pyplot.ioff()
job_server = pp.Server()

start = int(sys.argv[2])
stop = int(sys.argv[3])+1
for i in xrange(start,stop):
    fname = path + '{:0>3}'.format(i) + '.hdf5'
    snap = snapshot.load(fname)
    redshift = snap.header.Redshift
    for suffix in ['-dens-xy.png','-dens-xz.png','-dens-yz.png']:
        wpath = write_dir + '{:0>3}'.format(i) + suffix
        view = suffix[-6:-4]
        x,y,z = visualize.density(snap, view, boxsize, 1., length_unit, 
                                  job_server, pps, hsml_factor)
        print 'Plotting',view+'...'
        fig = pyplot.figure(1,(10,10))
        fig.clf()
        pyplot.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()])
        #pyplot.colorbar()
        ax = pyplot.gca()
        for sink in snap.sinks:
            if view == 'xy':
                ax.plot(sink[0],sink[1],'ko')
            elif view == 'xz':
                ax.plot(sink[0],sink[2],'ko')
            elif view == 'yz':
                ax.plot(sink[1],sink[2],'ko')
            else:
                raise KeyError
        ax.set_xlabel('AU')
        ax.set_ylabel('AU')
        ax.set_title('Redshift: %.5f' %redshift)
        pyplot.draw()
        pyplot.savefig(wpath, 
                       bbox_inches='tight')
    snap.close()
job_server.destroy()


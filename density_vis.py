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

length_unit = units.Length_AU
pps = 1000 # 'pixels' per side
hsml_factor = 1.7
path = os.getenv('HOME')+'/sim/vanilla/snapshot_'
write_dir = os.getenv('HOME')+'/data/simplots/vanilla/'
suffix = '-dens-xy.png'
view = suffix[-6:-4]
boxsize = 5e3
#pyplot.ioff()

job_server = pp.Server()
for i in [600]:
#for i in xrange(468,614):
    fname = path + '{:0>3}'.format(i) + '.hdf5'
    wpath = write_dir + '{:0>3}'.format(i) + suffix

    snap = snapshot.load(fname)
    redshift = snap.header.Redshift
    x,y,z = visualize.density(snap, view, boxsize, 1., length_unit, 
                              job_server, pps, hsml_factor)
    snap.close()
    print 'Plotting...'
    fig = pyplot.figure(1,(10,10))
    fig.clf()
    pyplot.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()])
    #pyplot.colorbar()
    ax = pyplot.gca()
    ax.set_xlabel('AU')
    ax.set_ylabel('AU')
    ax.set_title('Redshift: %.5f' %redshift)
    pyplot.draw()
    #pyplot.savefig(wpath, 
    #               bbox_inches='tight')
    #pyplot.show()

#!/usr/bin/env python
# gas_temp.py
# Jacob Hummel

import os
import sys
import glob
import numpy
from matplotlib import pyplot
#===============================================================================
if __name__ == '__main__':
    if ((len(sys.argv) not in [3,5]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gpe.py [simulation name] [snapshot]'
        print '   Option 2: python gpe.py [simulation 1] [snapshot 1] '\
              '[simulation 2] [snapshot 2]'
        sys.exit()

    write_dir = os.getenv('HOME')+'/data/simData/xrays/'
    simulation1 = sys.argv[1]
    snap1 = '{:0>4}'.format(sys.argv[2])
    fname = (os.getenv('HOME')+'/data/'+
             simulation1+'/haloz/'+snap1+'.npy')
    hdata1 = numpy.load(fname)
    if len(sys.argv) == 5:
        simulation2 = sys.argv[3]
        snap2 = '{:0>4}'.format(sys.argv[4])
        fname = (os.getenv('HOME')+'/data/'+
                 simulation2+'/haloz/'+snap2+'.npy')
        hdata2 = numpy.load(fname)

    fig = pyplot.figure(1,figsize=(10,10))
    fig.clf()
    ax = fig.gca()
    ax.set_xlabel(r'Radius [pc]')
    ax.set_ylabel(r'Gravitational Potential Energy [ergs]')
    sim1 = simulation1.split('/')[-1]
    ax.loglog(hdata1[:,1], hdata1[:,11], 'b-', label=sim1)
    if len(sys.argv) == 5:
        sim2 = simulation2.split('/')[-1]
        ax.loglog(hdata2[:,1], hdata2[:,11], 'r-', label=sim2)

    pyplot.legend(loc=2)
    ax.set_xlim(.1,200)
    #ax.set_ylim(5,1e6)
    if len(sys.argv) == 5:
        pyplot.savefig(write_dir+'gpe-'+sim1+'-'+sim2+'.png', 
                       bbox_inches='tight')
    else:
        pyplot.savefig(write_dir+'gpe-'+sim1+'.png', bbox_inches='tight')

#!/usr/bin/env python
# gas_temp.py
# Jacob Hummel

import os
import sys
import glob
import numpy
from matplotlib import pyplot,cm
from mpl_toolkits.mplot3d import Axes3D
#===============================================================================
if __name__ == '__main__':
    if ((len(sys.argv) not in [3]) or (sys.argv[1] == '-h')):
        print 'Usage :: python halo3d.py [property to plot] [simulation name]'
        sys.exit()
    keys = {'delta':(2, 'log$_{10}$ delta', 30, 40),
            'mass':(3, 'log$_{10}$ Mass [M$_{\odot}$]', 10, -50),
            'density':(4, 'log$_{10}$ Density [g cm$^{-3}$]'),
            'energy':(5, 'log$_{10}$ Energy [ergs]'),
            'Tshell':(6, 'log$_{10}$ Temperature [K]'),
            'Tavg':(7, 'log$_{10}$ Temperature [K]'),
            'nparticles':(8,'log$_{10}$ Number of Particles')}
    key = keys[sys.argv[1]]


    simulation = sys.argv[2]
    fname = (os.getenv('HOME')+'/data/'+
             simulation+'/halo_properties.npy')
    write_dir = os.getenv('HOME')+'/data/'+simulation+'/'
    hdata = numpy.load(fname)

    fig = pyplot.figure(1,figsize=(15,10))
    fig.clf()
    ax = fig.gca(projection='3d')
    hdata[hdata == 0.0] = numpy.nan
    hdata[:,:,key[0]] = numpy.log10(hdata[:,:,key[0]])
    haloz = numpy.ma.array(hdata, mask=numpy.isnan(hdata))
    ax.plot_surface(haloz[:,:,0],haloz[:,:,1],haloz[:,:,key[0]],
                    linewidth=.25, cmap=cm.autumn,cstride=3)

    ax.set_xlabel('Redshift')
    ax.set_ylabel('Radius [pc]')
    ax.set_zlabel(key[1])
    try:
        ax.view_init(elev=key[2],azim=key[3])
    except IndexError:
        pass
    pyplot.savefig(write_dir+'halo3d-'+sys.argv[1]+'.png', bbox_inches='tight')

#!/usr/bin/env python
# gas_phase.py
# Jacob Hummel

import os
import sys
import numpy
from matplotlib import pyplot
import pyGadget
#===============================================================================

def plot_phase(path, snap, write_dir, stride=50):
    wpath = write_dir + '{:0>4}'.format(snap)+'-phase.png'

    snapshot = pyGadget.snapshot.File(path)
    # Read relevant attributes
    h = snapshot.header.HubbleParam # H = 100*h
    redshift = snapshot.header.Redshift
    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    gamma = snapshot.gas.get_gamma()
    temp = snapshot.gas.get_temperature()
    h2frac = snapshot.gas.get_H2_fraction()
    electronfrac = snapshot.gas.get_electron_fraction()
    snapshot.close()
    
    # Initialization Complete --- Begin Analysis
    print 'Analyzing...'
    minimum = numpy.amin(particle_mass)
    refined = numpy.where(particle_mass <= minimum)[0][0:-1:stride]
    dens = dens[refined]
    gamma = gamma[refined]
    temp = temp[refined]
    h2frac = h2frac[refined]
    electronfrac = electronfrac[refined]

    ### Create Plot!
    ### All plots vs density
    print 'Plotting...'
    fig = pyplot.figure(1,(12,10))
    fig.clf()
    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)
    density_labels = (1e-2,1e0,1e2,1e4,1e6,1e8,1e10,1e12)

    # Temperature
    ax0.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xlim(1e-2, 1e12)
    ax0.set_ylim(7, 1e4)
    ax0.set_xticks(density_labels)
    ax0.set_xlabel('n [cm$^{-3}$]')
    ax0.set_ylabel('Temperature [K]')

    # Free electron fraction
    ax1.hexbin(dens,electronfrac,gridsize=500,bins='log',xscale='log',
               yscale='log',mincnt=1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1e-2, 1e12)
    ax1.set_ylim(1e-12, 1e-2)
    ax1.set_xticks(density_labels)
    ax1.set_xlabel('n [cm$^{-3}$]')
    ax1.set_ylabel('f$_{e^-}$')

    # Molecular Hydrogen Fraction
    ax2.hexbin(dens, h2frac,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1e-2, 1e12)
    ax2.set_ylim(1e-7,2)
    ax2.set_xticks(density_labels)
    ax2.set_xlabel('n [cm$^{-3}$]')
    ax2.set_ylabel('f$_{H_2}$')

    # Adiabatic exponent
    ax3.hexbin(dens, gamma,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax3.set_xscale('log')
    ax3.set_xlim(1e-2, 1e12)
    ax3.set_ylim(1,2)
    ax3.set_xticks(density_labels)
    ax3.set_xlabel('n [cm$^{-3}$]')
    ax3.set_ylabel('$\gamma_{ad}$')

    title = fig.suptitle('Redshift: %.3f' %(redshift,))
    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    pyplot.savefig(wpath, 
                   bbox_extra_artists=(title,),
                   bbox_inches='tight')
    return dens,temp
#===============================================================================

if __name__ == '__main__':

    pyplot.ioff()
    if len(sys.argv) < 4:
        print 'Usage: python gas_phase.py (simulation name) '\
            '(beginning snapshot) (final snapshot)'
        sys.exit()

    simulation = sys.argv[1]
    path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
    write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'

    start = int(sys.argv[2])
    stop = int(sys.argv[3])+1
    for snap in xrange(start,stop):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        dens, temp = plot_phase(fname, snap, write_dir, stride=50)

# particles3D.py
# plot particle positions in 3D
# Jacob Hummel

import os
import sys

import numpy
import h5py
from matplotlib import pyplot
import mpl_toolkits.axes_grid1 as axes_grid

from gadgetHDF5 import *
#===============================================================================

def plot_phase(path, write_dir, stride=50):
    snap = path[-8:-5]
    f = h5py.File(path,'r')
    
    # Code unit conversions
    unitMass_g = 1.989e43 # g
    unitLength_cm = 3.085678e21 # cm
    unitVelocity_cgs= 1.0e5
    unitDensity_cgs= unitMass_g / unitLength_cm**3
    unitTime_s= unitLength_cm / unitVelocity_cgs
    unitDensity_cgs= unitMass_g / unitLength_cm**3
    unitPressure_cgs= unitMass_g / unitLength_cm/ unitTime_s**2
    unitEnergy_cgs= unitMass_g * unitLength_cm**2 / unitTime_s**2
    # Fundamental constants
    k_B = 1.3806e-16 # erg/K
    m_H = 1.6726e-24 # g
    GRAVITY = 6.6726e-8 # dyne * cm**2 / g**2
    G = GRAVITY / unitLength_cm**3 * unitMass_g * unitTime_s**2
    X_h = 0.76 # Hydrogen Mass Fraction

    # Set hdf5 pathways
    header = f['Header']
    # Read relevant attributes
    boxSize = header.attrs.get('BoxSize')
    h = header.attrs.get("HubbleParam") # H = 100*h
    redshift = header.attrs.get('Redshift')
    f.close()

    particle_mass = partTypeN_readHDF5(path, 'Masses', 'gas')
    dens = partTypeN_readHDF5(path, 'Density', 'gas')
    energy =  partTypeN_readHDF5(path, 'InternalEnergy', 'gas')
    gamma =  partTypeN_readHDF5(path, 'Adiabatic index', 'gas')
    abundances = partTypeN_readHDF5(path, 'ChemicalAbundances', 'gas')

    # Initialization Complete --- Begin Analysis
    print 'Analyzing...'
    minimum = numpy.amin(particle_mass)
    refined = numpy.where(particle_mass <= minimum)[0][0:-1:stride]

    dens = dens[refined]
    energy = energy[refined]
    gamma = gamma[refined]
    abundances = abundances[refined]

    # Unit Conversions
    dens = dens * unitDensity_cgs * h*h * (1 + redshift)**3
    dens = dens * X_h / m_H
    energy = energy * unitEnergy_cgs / unitMass_g

    # Chemical Abundances
    # 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
    H2I = abundances[:,0]
    h2frac = 2*H2I
    mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
    mu = 1 / mu # mean molecular weight

    # Derived Properties
    temp = (mu * m_H / k_B) * energy * (gamma-1)
    #hot = numpy.where(temp > 8e2*dens**.5)[0]
    electronfrac = abundances[:,1] + abundances[:,4] + abundances[:,5]

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
    ax0.scatter(dens, temp, s=1, c='k', linewidths=0.0)
    #ax0.scatter(dens[hot], temp[hot], s=2, c='r', linewidths=0.0)
    ax0.set_xscale('log')
    ax0.set_yscale('log')

    ax0.set_xlim(1e-2, 1e12)
    ax0.set_ylim(7, 1e4)
    ax0.set_xticks(density_labels)
    ax0.set_xlabel('n [cm$^{-3}$]')
    ax0.set_ylabel('Temperature [K]')

    # Free electron fraction
    ax1.scatter(dens, electronfrac, s=1, c='k', linewidths=0.0)
    #ax1.scatter(dens[hot], electronfrac[hot], s=2, c='r', linewidths=0.0)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlim(1e-2, 1e12)
    ax1.set_ylim(1e-12, 1e-2)
    ax1.set_xticks(density_labels)
    ax1.set_xlabel('n [cm$^{-3}$]')
    ax1.set_ylabel('f$_{e^-}$')

    # Molecular Hydrogen Fraction
    ax2.scatter(dens, h2frac, s=1, c='k', linewidths=0.0)
    #ax2.scatter(dens[hot], h2frac[hot], s=2, c='r', linewidths=0.0)
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax2.set_xlim(1e-2, 1e12)
    ax2.set_ylim(1e-7,1)
    ax2.set_xticks(density_labels)
    ax2.set_xlabel('n [cm$^{-3}$]')
    ax2.set_ylabel('f$_{H_2}$')

    # Adiabatic exponent
    ax3.scatter(dens, gamma, s=1, c='k', linewidths=0.0)
    ax3.set_xscale('log')

    ax3.set_xlim(1e-2, 1e12)
    ax3.set_ylim(1,2)
    ax3.set_xticks(density_labels)
    ax3.set_xlabel('n [cm$^{-3}$]')
    ax3.set_ylabel('$\gamma_{ad}$')

    title = fig.suptitle('Redshift: %.3f' %(redshift,))
    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    pyplot.savefig(write_dir+snap+'-phase.png', 
                   bbox_extra_artists=(title,),
                   bbox_inches='tight')
#===============================================================================

if __name__ == '__main__':
    pyplot.ioff()
    wdir = os.getenv('HOME')+'/data/simplots/vanilla/'
    for snap in range(224,227): 
        path = (os.getenv('HOME')+'/sim/vanilla/snapshot_'+
                '{:0>3}'.format(snap)+'.hdf5')
        plot_phase(path, wdir, stride=50)

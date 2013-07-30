# particles3D.py
# plot particle positions in 3D
# Jacob Hummel

import os
import sys

import numpy
import h5py
from matplotlib import pyplot

from gadgetHDF5 import *
#===============================================================================

def plot_tvr(ax, path, stride=50):
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
    
    pos = partTypeN_readHDF5(path, 'Coordinates', 'gas')
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
    pos = pos[refined]
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    center = pos[dens.argmax()]

    # Unit Conversions
    dens = dens * unitDensity_cgs * h*h * (1 + redshift)**3
    dens = dens * X_h / m_H
    energy = energy * unitEnergy_cgs / unitMass_g
    pos = pos * unitLength_cm / 3.08e18 # convert to parsecs

    # Chemical Abundances
    # 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
    H2I = abundances[:,0]
    h2frac = 2*H2I
    mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
    mu = 1 / mu # mean molecular weight

    # Derived Properties
    temp = (mu * m_H / k_B) * energy * (gamma-1)
    radius = numpy.sqrt((x - center[0])**2 
                        + (y - center[1])**2
                        + (z - center[2])**2)
    hot = numpy.where(temp > 8e2*dens**.5)[0]

    print 'Plotting...'
    ax.scatter(radius, temp, s=.75, c='k', linewidths=0.0)
    ax.scatter(radius[hot], temp[hot], s=5, c='r', linewidths=0.0)
    return ax,redshift
    
#===============================================================================
def uniplot(path, write_dir, stride=50):
    snap = path[-8:-5]
    #Create Plot!
    fig = pyplot.figure(1,(15,10))
    fig.clf()
    ax = fig.add_subplot(111)
    ax,redshift = plot_tvr(ax, path,stride=stride)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(5e-2, 18)
    ax.set_ylim(10, 1e4)

    pyplot.title('Redshift: %.2f' %(redshift,))
    ax.set_xlabel('Radius [pc]')
    ax.set_ylabel('Temperature [K]')

    pyplot.savefig(write_dir+snap+'-tvr.png',
                   bbox_inches='tight')
#===============================================================================

if __name__ == '__main__':
    pyplot.ioff()
    wdir = os.getenv('HOME')+'/data/simplots/adaptive_gravsoft/'
    for snap in range(0,91): 
        path = (os.getenv('HOME')+'/sim/adaptive_gravsoft/snapshot_'+
                '{:0>3}'.format(snap)+'.hdf5')
        uniplot(path, wdir)

# particles3D.py
# plot particle positions in 3D
# Jacob Hummel

import os
import sys
import glob

import numpy
import h5py
import pp
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

from gadgetHDF5 import *
import zoomed
import directory

#===============================================================================
def plotgas(ax, path, mmh_com=(50,50,50), stride=1):
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
    a = 1/(1+redshift)
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
    refined = numpy.where(particle_mass <= minimum*8.0)[0][0:-1:stride]
    pos = pos[refined]
    particle_mass = particle_mass[refined]
    dens = dens[refined]
    energy = energy[refined]
    gamma = gamma[refined]
    abundances = abundances[refined]

    # Unit Conversions
    boxSize = boxSize/h*a
    pos = pos/h*a
    mmh_com = mmh_com/h*a

    dens = dens * unitDensity_cgs * h*h * (1 + redshift)**3
    dens = dens * X_h / m_H
    print 'Max Density: %.5e' %dens.max()
    energy = energy * unitEnergy_cgs / unitMass_g

    # 3D positions
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    x = x - mmh_com[0]
    y = y - mmh_com[1]
    z = z - mmh_com[2]

    L = boxSize/2
    x = numpy.where(x>L,(x-boxSize),x)
    x = numpy.where(x<-L,(x+boxSize),x)
    y = numpy.where(y>L,(y-boxSize),y)
    y = numpy.where(y<-L,(y+boxSize),y)
    z = numpy.where(z>L,(z-boxSize),z)
    z = numpy.where(z<-L,(z+boxSize),z)

    # Change units from kpc to pc
    x = x*1e3
    y = y*1e3
    z = z*1e3

    # Split in to most refined level and less refined
    refined = numpy.where(dens >= 7e9)[0]
    print refined.size, 'particles found'
    particle_mass = particle_mass[refined]
    dens = dens[refined]
    energy = energy[refined]
    gamma = gamma[refined]
    abundances = abundances[refined]
    x = x[refined]
    y = y[refined]
    z = z[refined]
    print ' x:: max: %.6e min: %.6e' %(x.max(),x.min())
    print ' y:: max: %.6e min: %.6e' %(y.max(),y.min())
    print ' z:: max: %.6e min: %.6e' %(z.max(),z.min())

    # Chemical Abundances
    # 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
    H2I = abundances[:,0]
    h2frac = 2*H2I
    mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
    mu = 1 / mu # mean molecular weight

    # Derived Properties
    temp = (mu * m_H / k_B) * energy * (gamma-1)
    hot = numpy.where(temp > 8e2*dens**.5)[0]

    print 'Plotting...'
    ax.scatter3D(x,y,z, s=10)
    return ax, redshift

#===============================================================================
def plotdm(ax, path, mmh_com=(50,50,50), stride=50):
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

    pos = partTypeN_readHDF5(path, 'Coordinates', 'dm')
    particle_mass = partTypeN_readHDF5(path, 'Masses', 'dm')

    # Initialization Complete --- Begin Analysis
    print 'Analyzing...'
    minimum = numpy.amin(particle_mass)
    refined = numpy.where(particle_mass <= minimum)[0][0:-1:stride]
    #unrefined = numpy.where(particle_mass <= minimum*8.0)[0][0:-1:stride]

    pos = pos[refined]
    particle_mass = particle_mass[refined]

    # Unit Conversions
    boxSize = boxSize/h
    pos = pos/h
    mmh_com = mmh_com/h

    # Initialization Complete --- Begin Analysis
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    x = x - mmh_com[0]
    y = y - mmh_com[1]
    z = z - mmh_com[2]

    L = boxSize/2
    x = numpy.where(x>L,(x-boxSize),x)
    x = numpy.where(x<-L,(x+boxSize),x)
    y = numpy.where(y>L,(y-boxSize),y)
    y = numpy.where(y<-L,(y+boxSize),y)
    z = numpy.where(z>L,(z-boxSize),z)
    z = numpy.where(z<-L,(z+boxSize),z)

    print 'Plotting...'
    ax.scatter3D(x, y, z, s=.5, c='k', linewidths=0.0)
    pyplot.draw()
    return ax, redshift

#===============================================================================
def uniplot(path, write_dir, ptype, stride=50):
    snap = path[-8:-5]
    #Create Plots!
    fig = pyplot.figure(1,(15,15))
    fig.clf()
    ax = fig.add_subplot(111,aspect='equal',projection='3d')
    if ptype == 'gas':
        ax,redshift = plotgas(ax, path, stride=stride)
    elif ptype == 'dm':
        ax,redshift = plotdm(ax, path, stride=stride)
    else:
        raise ValueError("ptype must be either 'dm' or 'gas'")
    ax.set_xlim3d(17.318,17.321)
    ax.set_ylim3d(13.9495,13.9525)
    ax.set_zlim3d(-20.305,-20.302)

    pyplot.title('Redshift: %.2f' %(redshift,))
    ax.set_xlabel('[physical pc]')
    ax.set_ylabel('[physical pc]')
    ax.set_zlabel('[physical pc]')
    #ax.mouse_init()
    pyplot.draw()
    pyplot.savefig(write_dir+snap+'-'+ptype+'.png',
                   bbox_inches='tight')

#===============================================================================

if __name__ == '__main__':
    #pyplot.ioff()
    wdir = os.getenv('HOME')+'/data/simplots/vanilla/'
    for snap in xrange(280,300): 
        path = (os.getenv('HOME')+'/sim/vanilla/snapshot_'+
                '{:0>3}'.format(snap)+'.hdf5')
        uniplot(path, wdir, 'gas', stride=1)
        #uniplot(path, wdir, 'dm', stride=50)

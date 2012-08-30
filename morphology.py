#!/usr/bin/env python
# morphology.py
# Jacob Hummel

import os
import sys
import numpy
from matplotlib import pyplot
from matplotlib.mlab import griddata

from mpl_toolkits.mplot3d import Axes3D
import pp

import pyGadget
#===============================================================================

length_unit = pyGadget.units.Length_AU
box = 1e3
depth = 200

pps = 1e3

write_dir = os.getenv('HOME')+'/data/simplots/vanilla-100/'
for snap in range(467,468):
    path = (os.getenv('HOME')+'/sim/vanilla-100/snapshot_'+
            '{:0>3}'.format(snap)+'.hdf5')

    snap = path[-8:-5]
    wpath = write_dir+snap+'-morphology.png'

    units = pyGadget.units
    constants = pyGadget.constants
    snapshot = pyGadget.snapshot.load(path)

    # Read relevant attributes
    h = snapshot.header.HubbleParam # H = 100*h
    redshift = snapshot.header.Redshift

    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    pos = snapshot.gas.get_coords(length_unit)
    snapshot.close()

    # Initialization Complete --- Begin Analysis
    print 'Analyzing...'
    # Select only highest resolution particles
    minimum = numpy.amin(particle_mass)
    refined = numpy.where(particle_mass <= minimum)[0]
    dens = dens[refined]
    pos = pos[refined]
    # Select only highest density particles
    #refined = numpy.where(dens < 1e8)[0]
    #dens = dens[refined]
    #pos = pos[refined]
    print 'Refinement complete.'
    
    center = pos[dens.argmax()]
    x = pos[:,0] - center[0]
    y = pos[:,1] - center[1]
    z = pos[:,2] - center[2]

    slice_ = numpy.where(numpy.abs(z) < depth/2)[0]
    dens = dens[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    slice_ = numpy.where(numpy.abs(x) < box/2)[0]
    dens = dens[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    slice_ = numpy.where(numpy.abs(y) < box/2)[0]
    dens = dens[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    print 'Slicing complete.'
    print ' x:: max: %.3e min: %.3e' %(x.max(),x.min())
    print ' y:: max: %.3e min: %.3e' %(y.max(),y.min())
    print ' density:: max: %.3e min: %.3e' %(dens.max(),dens.min())
    print ' Array size:', dens.size

    xres = yres = box/pps
    xvals = numpy.arange(-box/2,box/2,xres)
    yvals = numpy.arange(-box/2,box/2,yres)
    xi,yi = numpy.meshgrid(xvals,yvals)
    print 'Mesh created.'

    print 'Building image...'
    #points = [(i,j) for i in range(xvals.size) for j in range(yvals.size)]
    #print points

    print 'Triangulating...'
    zi = griddata(x,y,dens**2,xi,yi)
    print 'Triangulation complete.'

    print 'Plotting...'
    fig = pyplot.figure(1,(12,12))
    fig.clf()
    pyplot.imshow(zi)
    #pyplot.contourf(xi,yi,zi)
    #pyplot.contour(xi,yi,zi)
    #ax = fig.add_subplot(111,aspect='equal',projection='3d')
    #ax.scatter3D(x, y, z, s=5, c='k', linewidths=0.0)
    #ax.scatter3D(xi, yi, 0, s=1, c='r', linewidths=0.0)
    #ax.set_xlim3d(-box/2,box/2)
    #ax.set_ylim3d(-box/2,box/2)
    #ax.set_zlim3d(-box/2,box/2)

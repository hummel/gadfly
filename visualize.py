#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy
from matplotlib import pyplot
from matplotlib.mlab import griddata

import pp

import pyGadget
#===============================================================================

length_unit = pyGadget.units.Length_kpc
pps = 1000 # 'pixels' per side
hsml_factor = 1.7

write_dir = os.getenv('HOME')+'/data/simplots/vanilla-100/'
for snap in range(467,468):
    path = (os.getenv('HOME')+'/sim/vanilla-100/snapshot_'+
            '{:0>3}'.format(snap)+'.hdf5')

    snap = path[-8:-5]
    wpath = write_dir+snap+'-dens.png'

    units = pyGadget.units
    constants = pyGadget.constants
    snapshot = pyGadget.snapshot.load(path)

    # Read relevant attributes
    h = snapshot.header.HubbleParam
    a = snapshot.header.ScaleFactor
    redshift = snapshot.header.Redshift

    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    pos = snapshot.gas.get_coords(length_unit,no_h=False, comoving=True)
    smL = snapshot.gas.get_smoothing_length()
    sinkval = snapshot.gas.get_sinks()
    snapshot.close() 

    # Initialization Complete --- Begin Analysis
    print 'Analyzing...'
    # Select only highest resolution particles
    minimum = numpy.amin(particle_mass)
    refined = numpy.where(particle_mass <= minimum)[0]
    dens = dens[refined]
    smL = smL[refined]
    sinkval = sinkval[refined]
    pos = pos[refined]
    # Select only highest density particles
    #refined = numpy.where(dens < 1e8)[0]
    #dens = dens[refined]
    #pos = pos[refined]
    print 'Refinement complete.'

    width= 10*h/a/5e5
    depth= width
    center = pos[dens.argmax()]
    x = pos[:,0] - center[0]
    y = pos[:,1] - center[1]
    z = pos[:,2] - center[2]

    slice_ = numpy.where(numpy.abs(z) < depth/2)[0]
    dens = dens[slice_]
    smL = smL[slice_]
    sinkval = sinkval[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    slice_ = numpy.where(numpy.abs(x) < width/2)[0]
    dens = dens[slice_]
    smL = smL[slice_]
    sinkval = sinkval[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    slice_ = numpy.where(numpy.abs(y) < width/2)[0]
    dens = dens[slice_]
    smL = smL[slice_]
    sinkval = sinkval[slice_]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    print 'Slicing complete.'
    print ' x:: max: %.3e min: %.3e' %(x.max(),x.min())
    print ' y:: max: %.3e min: %.3e' %(y.max(),y.min())
    print ' density:: max: %.3e min: %.3e' %(dens.max(),dens.min())
    print ' sink values:: max: %.3e min: %.3e' %(sinkval.max(),sinkval.min())
    print ' smoothing length:: max: %.3e min: %.3e' %(smL.max(),smL.min())
    print ' Array size:', dens.size

    xres = yres = width/pps
    xvals = numpy.arange(-width/2,width/2,xres)
    yvals = numpy.arange(-width/2,width/2,yres)
    xi,yi = numpy.meshgrid(xvals,yvals)
    zi = numpy.zeros_like(xi)
    nzi = numpy.zeros_like(zi)



#===============================================================================
    zmin = 1e-2
    zmax = 1e1
    hsml = numpy.fmax(hsml_factor * smL, width / pps / 2.0)
    i_min = (x - hsml + width/2.0) / width*pps
    i_max = (x + hsml + width/2.0) / width*pps
    j_min = (y - hsml + width/2.0) / width*pps
    j_max = (y + hsml + width/2.0) / width*pps
    weight = dens*dens

    print 'Filling mesh...'
    for n in range(dens.size):
        percent = 100*n/dens.size
        if n%(dens.size/100)==0:
            print '%2i%%' %percent
        # if sink smoothing link is too inflated, 
        # artificially set it to accretion radius value
        if(sinkval[n] > 0):
            print 'sinkval > 0 !!!'
            hsml[n] = hsml_factor * 3.57101e-07

        for i in range(int(pps)):
            if(i >= i_min[n] and i <= i_max[n]):
                center_i = -width/2.0 + (i+0.5) * width/pps
                for j in range(int(pps)):
                    if(j >= j_min[n] and j <= j_max[n]):
                        center_j = -width/2.0 + (j+0.5) * width/pps
                        r2 = ((x[n] - center_i)**2
                              + (y[n] - center_j)**2) / hsml[n] / hsml[n]
                        if(r2 <= 1.0):
                            r = numpy.sqrt(r2)
                            if(r <= 0.5):
                                W_x = 1.0 - 6.0 * r**2 + 6.0 * r**3
                            else:
                                W_x = 2.0 * (1.0 - r)**3
                            zi[i][j] += weight[n] * dens[n] * W_x
                            nzi[i][j] += weight[n] * W_x

    zi = numpy.where(nzi > 0, zi/nzi,zi)
    zi = numpy.fmax(zi, zmin)
    zi = numpy.fmin(zi, zmax)
    zi = numpy.log10(zi)
    zi[0,0] = numpy.log10(zmin)
    zi[-1,-1] = numpy.log10(zmax)

#===============================================================================

              
    print 'Plotting...'
    fig = pyplot.figure(1,(10,10))
    fig.clf()
    pyplot.imshow(zi)
    #pyplot.contourf(xi,yi,zi)
    #pyplot.contour(xi,yi,zi)
    #ax = fig.add_subplot(111,aspect='equal',projection='3d')
    #ax.scatter3D(x, y, z, s=5, c='k', linewidths=0.0)
    #ax.scatter3D(xi, yi, 0, s=1, c='r', linewidths=0.0)
    #ax.set_xlim3d(-width/2,width/2)
    #ax.set_ylim3d(-width/2,width/2)
    #ax.set_zlim3d(-width/2,width/2)

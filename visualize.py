#!/usr/bin/env python
# visualize.py
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
pps = 1e3 # 'pixels' per side
hsml_factor = 1.7

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
    h = snapshot.header.HubbleParam
    a = snapshot.header.ScaleFactor
    redshift = snapshot.header.Redshift

    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    pos = snapshot.gas.get_coords(length_unit)
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

    width= 1e3
    depth= 200
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
    print 'Mesh created.'



#===============================================================================
    IDmass = 0.0
    IDsink = 0
    flag_i = 0
    flag_j = 0
    h = numpy.fmax(hsml_factor * smL, width / pps / 2.0)
    for n in range(dens.size):
        # if sink smoothing link is too inflated, 
        # artificially set it to accretion radius value
        if(sinkval[n] > 0):
            h[n] = hsml_factor * 3.57101e-07

        weight = dens*dens
        i = 0
        j=0
        i_min = int((x[n] - h[n] + width / 2.0) / width * pps)
        if(i_min < 0):
            i_min = 0
        i_max = int((x[n] + h[n] + width / 2.0) / width * pps)
        if(i_max > pps-1):
            i_max = pps-1
        j_min = int((y[n] - h[n] + width / 2.0) / width * pps)
        if(j_min < 0):
            j_min = 0
        j_max = int((y[n] + h[n] + width / 2.0) / width * pps)
        if(j_max > pps-1):
            j_max = pps-1
        print i_min, i_max, j_min, j_max
    """
	      do
		
		  if(i >= i_min && i <= i_max)
		    
		      flag_i = 1
		      
		      center_i = center_x - width / 2.0 + (i + 0.5) * width / (double) pps
		      do
			
			  if(j >= j_min && j <= j_max)
			    
			      flag_j = 1
			      
			      center_j = center_y - width / 2.0 + (j + 0.5) * width / (double) pps
			      
			      x2 = ((P[n].Pos[0] - center_i) * (P[n].Pos[0] - center_i) + (P[n].Pos[1] - center_j) * (P[n].Pos[1] - center_j)) / h / h
			      
			      if(x2 <= 1.0)
				
				  x = sqrt(x2)
				  
				  if(x <= 0.5)
				    W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x
				  else
				    W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x)
				  
				  grid1[i][j] += weight * P[n].nh * W_x
				  
				  n_grid1[i][j] += weight * W_x
				
			      
			      
			      
			    
			  else if(j > j_max)
			    
			      flag_j = 2
			    
			  else
			    
			      flag_j = 0
			    
			  
			  j++
			
		      while((flag_j == 0 || flag_j == 1) && j < pps)
		      
		      j = 0
		    
		  else if(i > i_max)
		    
		      flag_i = 2
		    
		  else
		    
		      flag_i = 0
		    
		  
		  i++
		
	      while((flag_i == 0 || flag_i == 1) && i < pps)
	      
	      i = 0
	    
	
      
      for(i = 0 i < pps i++)
	
	  for(j = 0 j < pps j++)
	    
	      if(n_grid1[i][j] > 0)
		
		  grid1[i][j] /= double(n_grid1[i][j])
		
	      
	      if(grid1[i][j] < min)
		
		  grid1[i][j] = min
		
	      
	      if(grid1[i][j] > max)
		
		  grid1[i][j] = max
		
	      
	      grid1[i][j] = log10(grid1[i][j])
	    
	

#===============================================================================

              
    """
    #print 'Plotting...'
    #fig = pyplot.figure(1,(12,12))
    #fig.clf()
    #pyplot.imshow(zi)
    #pyplot.contourf(xi,yi,zi)
    #pyplot.contour(xi,yi,zi)
    #ax = fig.add_subplot(111,aspect='equal',projection='3d')
    #ax.scatter3D(x, y, z, s=5, c='k', linewidths=0.0)
    #ax.scatter3D(xi, yi, 0, s=1, c='r', linewidths=0.0)
    #ax.set_xlim3d(-width/2,width/2)
    #ax.set_ylim3d(-width/2,width/2)
    #ax.set_zlim3d(-width/2,width/2)

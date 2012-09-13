#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy

from scipy import weave
from matplotlib import pyplot
from matplotlib.mlab import griddata
import pp

import pyGadget
#===============================================================================
def scalar_map(pps,width, x,y,scalar_field,hsml,zshape):
    zi = numpy.zeros(zshape)
    nzi = numpy.zeros_like(zi)

#/##==== C code ====
    code = \
        """
        # line 23 visualize.py
        int i,j, i_min,i_max,j_min,j_max;
        double r,r2,weight,W_x;
        for(n =0; n < N_gas; n++) 
          {
            i = 0;
            j = 0;
            i_min = int((x[n] - hsml[n] + width/2.0) / width*pps);
            i_max = int((x[n] + hsml[n] + width/2.0) / width*pps);
            j_min = int((y[n] - hsml[n] + width/2.0) / width*pps);
            j_max = int((y[n] + hsml[n] + width/2.0) / width*pps);
            weight = scalar_field*scalar_field;
            do 
              {
                if(i >= i_min && i <= i_max[n])
                  {
                    flag_i = 1;
                    center_i = -width/2.0 + (i+0.5) * width/ (double) pps;
                    do
                      {
                        if(j >= j_min && j <= j_max)
                          {
                            flag_j = 1;
                            center_j = -width/2.0 + (j+0.5)
                              * width / (double) pps;
                            r2 = ((x[n] - center_i) * (x[n] - center_i)
                                  + (y[n] - center_j) * (y[n] - center_j))
                              / hsml[n] / hsml[n];
                            if(r2 <= 1.0)
                              {
                                r = sqrt(r2);
                                if(r <= 0.5)
                                  W_x = 1.0 - 6.0 * r*r + 6.0 * r*r*r;
                                else
                                  W_x = 2.0 * (1.0-r) * (1.0-r) * (1.0-r);
                                zi[i][j] += weight[n] * scalar_field[n] * W_x;
                                nzi[i][j] += weight[n] * W_x;
                              }
                          }
                        else if(j > j_max)
                          {
                            flag_j = 2;
                          }
                        else
                          {
                            flag_j = 0;
                          }
                        j++;
                      }
                    while((flag_j == 0 || flag_j == 1) && j < N_grid);
                    j = 0;
                  }
                else if(i > i_max)
                  {
                    flag_i = 2;
                  }
                else
                  {
                    flag_i = 0;
                  }
                i++;
              }
            while((flag_i == 0 || flag_i == 1) && i < N_grid);
            i = 0;
          }
        return_val = zi
        """
    zi = weave.inline(code,['pps','width',' x','y',
                            'scalar_field','hsml','zi','nzi'],
                      verbose=2,force=1)
    return zi,nzi
#===============================================================================

length_unit = pyGadget.units.Length_kpc
pps = 800 # 'pixels' per side
hsml_factor = 1.7

write_dir = os.getenv('HOME')+'/data/simplots/vanilla-100/'
for snap in xrange(467,468):
    #Create jobserver
    job_server = pp.Server()#ppservers=("*",))
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

    width= 1e-4
    depth= width/10
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
    zshape = zi.shape

    hsml = numpy.fmax(hsml_factor * smL, width / pps / 2.0)
    # if sink smoothing link is too inflated, 
    # artificially set it to accretion radius value
    #if(sinkval[n] > 0):
    #    print 'sinkval > 0 !!!'
    #    hsml[n] = hsml_factor * 3.57101e-07
    """
    print 'Distributing...'
    jobs = []
    server_list = job_server.get_active_nodes()
    ncpus = sum(server_list.values())
    #Divide in to 8X as many tasks as there are cpus.
    parts = ncpus*8
    start = 0
    end = dens.size - 1
    step = (end - start) / parts + 1
    print 'Dividing the input into',parts,'arrays of',step,'particles'
    print 'for calculation on',ncpus,'cpus.'
    for cpu in xrange(parts):
        nstart = start+cpu*step
        nend = min(start+(cpu+1)*step, end)
        #'''
        jobs.append(job_server.submit(scalar_map,
                                      (pps,width, 
                                       x[nstart:nend],
                                       y[nstart:nend],
                                       dens[nstart:nend],
                                       hsml[nstart:nend],zshape),
                                      (),('numpy','from scipy import weave')))
        #'''
    #sys.exit()
    print 'Calculating...'
    for job in jobs:
        pzi,pnzi = job()
        zi += pzi
        nzi += pnzi

    job_server.print_stats()
    job_server.destroy()
    """
    zi,nzi = scalar_map(pps,width,x,y,dens,hsml,zshape),
    zi = numpy.where(nzi > 0, zi/nzi, zi)
    #zi = numpy.fmax(zi, zmin)
    #zi = numpy.fmin(zi, zmax)
    zi = numpy.log10(zi)
    #zi[0,0] = numpy.log10(zmin)
    #zi[-1,-1] = numpy.log10(zmax)

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
    pyplot.show()

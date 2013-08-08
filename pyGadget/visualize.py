#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy

from scipy import weave
from scipy.weave import converters

from . import analyze
#===============================================================================
def scalar_map(pps,width, x,y,scalar_field,hsml,zshape):
    zi = numpy.zeros(zshape)
    nzi = numpy.zeros_like(zi)
    N_gas = scalar_field.size
    code = \
        r"""
        int procs, nthreads;
        int i,j,n, i_min,i_max,j_min,j_max;
        int flag_i = 0;
        int flag_j = 0;
        double center_i,center_j;
        double r,r2,weight,W_x;

        /* Get environment information */
        procs = omp_get_num_procs();
        nthreads = omp_get_num_threads();
        /* Print environment information */
        printf("Number of processors = %d\n", procs);
                
        #pragma omp parallel for \
          private(n,i,j,i_min,i_max,j_min,j_max,flag_i,flag_j, \
          center_i,center_j,r,r2,weight,W_x)
        for(n =0; n < N_gas; n++) 
          {
            i = 0;
            j = 0;
            i_min = int((x(n) - hsml(n) + width/2.0) / width*pps);
            i_max = int((x(n) + hsml(n) + width/2.0) / width*pps);
            j_min = int((y(n) - hsml(n) + width/2.0) / width*pps);
            j_max = int((y(n) + hsml(n) + width/2.0) / width*pps);
            weight = scalar_field(n)*scalar_field(n);
            do 
              {
                if(i >= i_min && i <= i_max)
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
                            r2 = ((x(n) - center_i) * (x(n) - center_i)
                                  + (y(n) - center_j) * (y(n) - center_j))
                              / hsml(n) / hsml(n);
                            if(r2 <= 1.0)
                              {
                                r = sqrt(r2);
                                if(r <= 0.5)
                                  W_x = 1.0 - 6.0 * r*r + 6.0 * r*r*r;
                                else
                                  W_x = 2.0 * (1.0-r) * (1.0-r) * (1.0-r);
                                zi(i,j) += weight * scalar_field(n) * W_x;
                                nzi(i,j) += weight * W_x;
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
                    while((flag_j == 0 || flag_j == 1) && j < pps);
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
            while((flag_i == 0 || flag_i == 1) && i < pps);
            i = 0;
          }
        """
    weave.inline(code,
                 ['pps','width','x','y','scalar_field',
                  'hsml','zi','nzi','N_gas'],
                 compiler='gcc',
                 headers=['<stdio.h>','<math.h>','<omp.h>'],
                 extra_compile_args=['-fopenmp ' ],
                 libraries=['gomp'], type_converters=converters.blitz)
    zi = numpy.where(nzi > 0, zi/nzi, zi)
    return zi

#===============================================================================
def py_scalar_map(pps,width, x,y,scalar_field,hsml,zshape):
    zi = numpy.zeros(zshape)
    nzi = numpy.zeros_like(zi)
    i_min = (x - hsml + width/2.0) / width*pps
    i_max = (x + hsml + width/2.0) / width*pps
    j_min = (y - hsml + width/2.0) / width*pps
    j_max = (y + hsml + width/2.0) / width*pps
    weight = scalar_field*scalar_field
    for n in range(scalar_field.size):
        for i in xrange(pps):
            if(i >= i_min[n] and i <= i_max[n]):
                center_i = -width/2.0 + (i+0.5) * width/pps
                for j in xrange(pps):
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
                            zi[i][j] += weight[n] * scalar_field[n] * W_x
                            nzi[i][j] += weight[n] * W_x
    zi = numpy.where(nzi > 0, zi/nzi, zi)
    return zi

#===============================================================================
def set_viewpoint(pos, dens, viewpoint, centering, **kwargs):
    if viewpoint == 'xy':
        x = pos[:,0]
        y = pos[:,1]
        z = pos[:,2]
    elif viewpoint == 'xz':
        x = pos[:,0]
        y = pos[:,2]
        z = pos[:,1]
    elif viewpoint == 'yz':
        x = pos[:,1]
        y = pos[:,2]
        z = pos[:,0]
    else:
        raise KeyError

    if centering == 'halo':
        x,y,z = analyze.find_center(x,y,z,dens,**kwargs)

    elif centering == 'box':
        x = x - (x.max() + x.min())/2
        y = y - (y.max() + y.min())/2
        z = z - (z.max() + z.min())/2
    else:
        raise KeyError
    return x,y,z

#===============================================================================
def density_projection(snapshot, width, thickness, length_unit, *args,**kwargs):
    pps = kwargs.pop('pps', 500)
    sm = kwargs.pop('sm', 1.7)
    # Read relevant attributes
    h = snapshot.header.HubbleParam
    a = snapshot.header.ScaleFactor
    redshift = snapshot.header.Redshift
    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    pos = snapshot.gas.get_coords(length_unit)
    smL = snapshot.gas.get_smoothing_length(length_unit)
    sinkval = snapshot.gas.get_sinks()

    x,y,z = set_viewpoint(pos,dens,*args,**kwargs)

    # Save sink particle positions for overplotting
    snapshot.sinks = []
    for sink_id in snapshot.sink_ids:
        snapshot.sinks.append((x[sink_id],y[sink_id],z[sink_id],
                               particle_mass[sink_id]))
    try: 
        assert dens.max() <= 1e12
    except AssertionError: 
        print 'Warning: Maximum density exceeds 1e12 particles/cc!'
        print 'Max Density: %.5e' %dens.max()
    if thickness:
        depth= width*thickness
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
    print ' x:: max: %.3e min: %.3e' %(x.max(),x.min())
    print ' y:: max: %.3e min: %.3e' %(y.max(),y.min())
    print ' density:: max: %.3e min: %.3e' %(dens.max(),dens.min())
    print ' Array size:', dens.size

    xres = yres = width/pps
    xvals = numpy.arange(-width/2,width/2,xres)
    yvals = numpy.arange(-width/2,width/2,yres)
    xi,yi = numpy.meshgrid(xvals,yvals)
    zi = numpy.zeros_like(xi)
    nzi = numpy.zeros_like(zi)
    zshape = zi.shape
    hsml = numpy.fmax(sm * smL, width / pps / 2.0)

    print 'Calculating...'
    zi = scalar_map(pps,width,x,y,dens,hsml,zshape)
    print 'density:: min: %.3e max: %.3e' %(zi.min(),zi.max())
    return xi,yi,zi
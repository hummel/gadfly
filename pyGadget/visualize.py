#!/usr/bin/env python
# visualize.py
# Jacob Hummel

import os
import sys
import numpy
from numba import jit
from scipy import weave
from scipy.weave import converters

import analyze
import coordinates
#===============================================================================
def scalar_map(y,x,scalar_field,hsml,width,pps,zshape):
    """
    code contains an algorithm for doing sph particle smoothing.
    """
    # !!!!! x and y are reversed in the c-code below - no idea why/how...
    # !!!!! (Hence the reversal above)
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
def py_scalar_map(y,x,scalar_field,hsml,width,pps,zshape):
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
@jit(nopython=True, nogil=True)
def numba_loop(y,x,scalar_field,hsml,zi,nzi,weight,npart,width,pps,i_min,i_max,j_min,j_max):
    for n in xrange(npart):
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
    return zi, nzi
def numba_scalar_map(y,x,scalar_field,hsml,width,pps,zshape):
    zi = numpy.zeros(zshape)
    nzi = numpy.zeros_like(zi)
    i_min = (x - hsml + width/2.0) / width*pps
    i_max = (x + hsml + width/2.0) / width*pps
    j_min = (y - hsml + width/2.0) / width*pps
    j_max = (y + hsml + width/2.0) / width*pps
    weight = scalar_field*scalar_field
    npart = scalar_field.size
    zi,nzi = numba_loop(y,x,hsml,scalar_field,zi,nzi,weight,
                        npart,width,pps,i_min,i_max,j_min,j_max)
    zi = numpy.where(nzi > 0, zi/nzi, zi)
    return zi

#===============================================================================
def set_view(view, xyz, **kwargs):
    uvw = kwargs.pop('velocity', None)
    dens = kwargs.pop('density', None)
    mass = kwargs.pop('mass', None)
    dlim = kwargs.pop('dens_lim', 1e9)
    if view ==  'face':
        if uvw is None or dens is None:
            raise KeyError("setting face-on view requires density, "\
                           "position and velocity!")
        pos, vel = analyze.data_slice(dens > dlim, xyz, uvw)
        axis, angle = analyze.faceon_rotation(pos, vel, mass)
        xyz = coordinates.rotate(xyz, axis, angle)
        uvw = coordinates.rotate(uvw, axis, angle)
    elif view == 'xy':
        pass
    elif view == 'xz':
        xyz = coordinates.rotate(xyz,'x', numpy.pi/2)
        if uvw is not None:
            uvw = coordinates.rotate(uvw,'x', numpy.pi/2)
    elif view == 'yz':
        xyz = coordinates.rotate(xyz,'z', numpy.pi/2)
        xyz = coordinates.rotate(xyz,'x', numpy.pi/2)
        if uvw is not None:
            uvw = coordinates.rotate(uvw,'z', numpy.pi/2)
            uvw = coordinates.rotate(uvw,'x', numpy.pi/2)
    else:
        for rot in view:
            xyz = coordinates.rotate(xyz, rot[0], rot[1])
        if uvw is not None:
            for rot in view:
                uvw = coordinates.rotate(uvw, rot[0], rot[1])
    if uvw is not None:
        return xyz, uvw
    else:
        return xyz

def trim_view(width, x, y, z, *args, **kwargs):
    depth = kwargs.pop('depth',1.0)
    arrs = [i for i in args]
    if depth:
        depth *= width
        slice_ = numpy.where(numpy.abs(z) < depth/2)[0]
        x = x[slice_]
        y = y[slice_]
        z = z[slice_]
        for i,arr in enumerate(arrs):
            arrs[i] = arr[slice_]
    slice_ = numpy.where(numpy.abs(x) < width/2)[0]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    for i,arr in enumerate(arrs):
        arrs[i] = arr[slice_]
    slice_ = numpy.where(numpy.abs(y) < width/2)[0]
    x = x[slice_]
    y = y[slice_]
    z = z[slice_]
    for i,arr in enumerate(arrs):
        arrs[i] = arr[slice_]
    print ' x:: max: %.3e min: %.3e' %(x.max(),x.min())
    print ' y:: max: %.3e min: %.3e' %(y.max(),y.min())
    return [x,y,z]+arrs

def build_grid(width,pps):
    xres = yres = width/pps
    xvals = numpy.arange(-width/2,width/2,xres)
    yvals = numpy.arange(-width/2,width/2,yres)
    return numpy.meshgrid(xvals,yvals)

def project(snapshot, loadable, scale, view, **kwargs):
    pps = kwargs.pop('pps',500)
    sm = kwargs.pop('sm',1.7)
    shiftx = kwargs.pop('shiftx',None)
    shifty = kwargs.pop('shifty',None)
    shiftz = kwargs.pop('shiftz',None)
    dens_lim = kwargs.pop('dens_lim', None)
    boxsize = float("".join(ch if ch.isdigit() else "" for ch in scale))
    unit = "".join(ch if not ch.isdigit() else "" for ch in scale)
    scalar = snapshot.gas._load_dict[loadable]()
    pos = snapshot.gas.get_coords(unit)
    hsml = snapshot.gas.get_smoothing_length(unit)
    if loadable not in ['density', 'ndensity']:
        dens = snapshot.gas.get_number_density()
    else:
        dens = scalar

    print 'Calculating...'
    pos = analyze.center_box(pos,density=dens,**kwargs)
    pos = set_view(view, pos)
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    if shiftx:
        x += shiftx
    if shifty:
        y += shifty
    if shiftz:
        z += shiftz
    snapshot.update_sink_coordinates(x,y,z)
    # Artificially shrink sink smoothing lengths.
    for s in snapshot.sinks:
        hsml[s.index] *= .5
    arrs = [scalar, hsml]
    if loadable in ['density', 'ndensity']:
        x,y,z,scalar,hsml = trim_view(boxsize, x,y,z,scalar,hsml,**kwargs)
        if dens_lim:
            arrs = [scalar,x,y,z,hsml]
            scalar,x,y,z,hsml = analyze.data_slice(scalar > dens_lim, *arrs)
    else:
        arrs.append(dens)
        x,y,z,scalar,hsml,dens = trim_view(boxsize, x,y,z, *arrs, **kwargs)
        if dens_lim:
            arrs = [dens,x,y,z,scalar,hsml]
            dens,x,y,z,scalar,hsml = analyze.data_slice(dens > dens_lim, *arrs)
    hsml = numpy.fmax(sm * hsml, boxsize/pps/2)
    xi,yi = build_grid(boxsize,pps)
    zi = scalar_map(x,y,scalar,hsml,boxsize,pps,xi.shape)
    print '%s:: min: %.3e max: %.3e' %(loadable, zi.min(),zi.max())
    imscale = kwargs.pop('imscale','log')
    if imscale == 'log':
        zi = numpy.log10(zi)
        print 'log(%s):: min: %.3e max: %.3e' %(loadable, zi.min(),zi.max())
    else:
        print 'Returning raw (non-log) values!'
    return xi,yi,zi

# analyze.py
# Jacob Hummel
import numpy
import units
import constants

#===============================================================================
def reject_outliers(data, m=2):
    return data[numpy.abs(data - numpy.mean(data)) < m * numpy.std(data)]

def find_center(x, y, z, dens=None, **kwargs):
    centering = kwargs.pop('centering','box')
    verbose = kwargs.get('verbose', True)
    if centering in ['avg', 'max']:
        if dens is not None:
            dens_limit = kwargs.pop('dens_limit', 1e11)
            nparticles = kwargs.pop('centering_npart', 10)
            if centering == 'avg':
                hidens = numpy.where(dens >= dens_limit)[0]
                while hidens.size < nparticles:
                    dens_limit /= 2
                    hidens = numpy.where(dens >= dens_limit)[0]
                if verbose:
                    print ('Center averaged over %d particles' %nparticles)
                    print ('Center averaged over all particles with density '\
                               'greater than %.2e particles/cc' %dens_limit)
                #Center on highest density clump, rejecting outliers:
                cx = numpy.average(reject_outliers(x[hidens]))
                cy = numpy.average(reject_outliers(y[hidens]))
                cz = numpy.average(reject_outliers(z[hidens]))
                print 'Density averaged box center: %.3e %.3e %.3e' %(cx,cy,cz)
            elif centering == 'max':
                center = dens.argmax()
                cx,cy,cz = x[center], y[center], z[center]
                print 'Density maximum box center: %.3e %.3e %.3e' %(cx,cy,cz)
        else:
            raise KeyError("'avg' and 'max' centering require gas density")
    elif centering == 'box':
        cx = (x.max() + x.min())/2
        cy = (y.max() + y.min())/2
        cz = (z.max() + z.min())/2
        print 'Simple box center: %.3e %.3e %.3e' %(cx,cy,cz)
    else:
        raise KeyError("Centering options are 'avg', 'max', and 'box'")
    return cx,cy,cz

def center_box(pos, center=None, **kwargs):
    dens = kwargs.pop('density', None)
    centering = kwargs.get('centering', None)
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    if center:
        cx = center[0]
        cy = center[1]
        cz = center[2]
    elif centering:
        cx,cy,cz = find_center(x,y,z, dens, **kwargs)
    else:
        print "WARNING! NO CENTER OR CENTERING ALGORITHM SPECIFIED!"
        print "Attempting simple box centering..."
        cx,cy,cz = find_center(x,y,z, centering='box')

    x -= cx
    y -= cy
    z -= cz
    return numpy.column_stack((x,y,z))

def cart2sph(x, y, z):
    r = numpy.sqrt(numpy.square(x) + numpy.square(y) + numpy.square(z))
    theta = numpy.arccos(z/r)
    phi = numpy.arctan2(y,x)
    return r,theta,phi

def sph2cart(r, theta, phi):
    x = r * numpy.sin(theta) * numpy.cos(phi)
    y = r * numpy.sin(theta) * numpy.sin(phi)
    z = r * numpy.cos(theta)
    return x,y,z

def cart2cyl(x ,y, z):
    r = numpy.sqrt(numpy.square(x) + numpy.square(y))
    theta = numpy.arctan2(y,x)
    return r,theta,z

def cyl2cart(r, theta, z):
    x = r * numpy.cos(theta)
    y = r * numpy.sin(theta)
    return x,y,z

def cartesian_unit_vectors(xyz):
    unit_x = numpy.zeros_like(xyz)
    unit_y = numpy.zeros_like(xyz)
    unit_z = numpy.zeros_like(xyz)
    unit_x[:,0] = 1.
    unit_y[:,1] = 1.
    unit_z[:,2] = 1.
    return unit_x, unit_y, unit_z

def spherical_unit_vectors(xyz):
    unit_x, unit_y, unit_z = cartesian_unit_vectors(xyz)
    r,theta,phi = cart2sph(xyz[:,0],xyz[:,1],xyz[:,2])
    r = r[:,numpy.newaxis]
    theta = theta[:,numpy.newaxis]
    phi = phi[:, numpy.newaxis]
    unit_r = xyz / r
    unit_phi = unit_y * numpy.cos(phi) - unit_x * numpy.sin(phi)
    unit_theta = (unit_x * numpy.cos(theta) * numpy.cos(phi)
                  + unit_y * numpy.cos(theta) * numpy.sin(phi)
                  - unit_z * numpy.sin(theta))
    return unit_r, unit_theta, unit_phi

def cart2sph_velocities(xyz, vxyz):
    unit_r, unit_theta, unit_phi = spherical_unit_vectors(xyz)
    # numpy.einsum('ij,ij->i',...) is einstein notation for dot product.
    vr = numpy.einsum('ij,ij->i',vxyz,unit_r)
    vtheta = numpy.einsum('ij,ij->i',vxyz,unit_theta)
    vphi = numpy.einsum('ij,ij->i',vxyz,unit_phi)
    return vr, vtheta, vphi

def density_cut(dens_lim, dens, *args):
    arrs = [i for i in args]
    slice_ = numpy.where(numpy.abs(dens) > dens_lim)[0]
    dens = dens[slice_]
    for i,array in enumerate(arrs):
        arrs[i] = array[slice_]
    print ' density cut:: max: %.3e min: %.3e' %(dens.max(),dens.min())
    return [dens]+arrs

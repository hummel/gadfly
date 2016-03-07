# analyze.py
# Jacob Hummel
import numpy
import pandas
import units

#===============================================================================
def reject_outliers(data, m=2):
    return data[numpy.abs(data - numpy.mean(data)) < m * numpy.std(data)]

def find_center(pos_vel, density=None, **kwargs):
    centering = kwargs.pop('centering','box')
    verbose = kwargs.get('verbose', True)
    if centering in ['avg', 'max']:
        if density is not None:
            dens_limit = kwargs.pop('dens_limit', 1e8)
            nparticles = kwargs.pop('centering_npart', 100)
            if centering == 'avg':
                hidens = numpy.where(density >= dens_limit)[0]
                while hidens.size < nparticles:
                    dens_limit /= 2
                    hidens = numpy.where(density >= dens_limit)[0]
                if verbose:
                    print ('Center averaged over %d particles' %hidens.size)
                    print ('Center averaged over all particles with density '\
                               'greater than %.2e particles/cc' %dens_limit)
                #Center on highest density clump, rejecting outliers:
                center = reject_outliers(pos_vel).mean()
                print 'Density averaged box center:',
            elif centering == 'max':
                center = pos_vel.iloc[density.argmax()]
                print 'Density maximum box center:',
        else:
            raise KeyError("'avg' and 'max' centering require gas density")
    elif centering == 'box':
        center = (pos_vel.max() + pos_vel.min())/2
        try:
            center[['u', 'v', 'w']] = 0
        except ValueError:
            pass
        print 'Simple box center:',
    else:
        raise KeyError("Centering options are 'avg', 'max', and 'box'")
    print '%.3e %.3e %.3e' %(center.x, center.y, center.z)
    return center

def center_box(pos_vel, center=None, vcenter=None, **kwargs):
    density = kwargs.pop('density', None)
    centering = kwargs.get('centering', None)
    if center:
        center = pandas.Series(center, index=['x', 'y', 'z'])
        vc = pandas.Series([0,0,0], index=['u', 'v', 'w'])
        center = pandas.concat([center, vc])
        if vcenter is not None:
            center['u'] = vcenter[0]
            center['v'] = vcenter[1]
            center['w'] = vcenter[2]
    elif centering:
        center = find_center(pos_vel, density, **kwargs)
    else:
        print "WARNING! NO CENTER OR CENTERING ALGORITHM SPECIFIED!"
        print "Attempting simple box centering..."
        center = find_center(pos_vel, density, **kwargs)
    pos_vel -= center
    return pos_vel

def angular_momentum(xyz, uvw, mass):
    rxv = numpy.cross(xyz, uvw)
    L = (mass[:, numpy.newaxis] * rxv)
    return L

def total_angular_momentum(xyz, uvw, mass, L=None):
    if L is None:
        L = angular_momentum(xyz, uvw, mass)
    return L.sum(axis=0)

def moment_of_inertia(xyz, uvw, mass, L=None):
    if L is None:
        L = total_angular_momentum(xyz, uvw, mass)
    elif L.size > 3:
        L = total_angular_momentum(xyz, uvw, mass, L)
    unitL = L / numpy.linalg.norm(L)
    rxL = numpy.cross(xyz, unitL)
    rxL2 = numpy.einsum('ij,ij->i',rxL,rxL)
    I = (mass[:, numpy.newaxis] * rxL2).sum()
    return I

def angular_velocity(xyz, uvw, mass, L=None, I=None):
    if L is None:
        L = total_angular_momentum(xyz, uvw, mass)
    if I is None:
        I = moment_of_inertia(xyz, uvw, mass, L)
    return L/I

def faceon_rotation(xyz, uvw, mass=None):
    if mass is None:
        mass = numpy.ones(xyz.shape[0])
    L = total_angular_momentum(xyz, uvw, mass)
    unitL = L / numpy.linalg.norm(L)
    z = numpy.array([0.,0.,1.])
    axis = numpy.cross(z, unitL)
    angle = numpy.arccos(unitL.dot(z))
    return axis, angle

# analyze.py
# Jacob Hummel
import numpy
import units
import constants

#===============================================================================
def reject_outliers(data, m=2):
    return data[numpy.abs(data - numpy.mean(data)) < m * numpy.std(data)]

def data_slice(expr, *args):
    arrs = [i for i in args]
    slice_ = numpy.where(expr)[0]
    for i,array in enumerate(arrs):
        arrs[i] = array[slice_]
    return arrs

def find_center(x, y, z, dens=None, **kwargs):
    centering = kwargs.pop('centering','box')
    verbose = kwargs.get('verbose', True)
    if centering in ['avg', 'max']:
        if dens is not None:
            dens_limit = kwargs.pop('dens_limit', 1e8)
            nparticles = kwargs.pop('centering_npart', 100)
            if centering == 'avg':
                hidens = numpy.where(dens >= dens_limit)[0]
                while hidens.size < nparticles:
                    dens_limit /= 2
                    hidens = numpy.where(dens >= dens_limit)[0]
                if verbose:
                    print ('Center averaged over %d particles' %hidens.size)
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
    return (cx, cy, cz)

def find_center_vcenter(x, y, z, u, v, w, dens=None, **kwargs):
    centering = kwargs.pop('centering','box')
    verbose = kwargs.get('verbose', True)
    if centering in ['avg', 'max']:
        if dens is not None:
            dens_limit = kwargs.pop('dens_limit', 1e8)
            nparticles = kwargs.pop('centering_npart', 100)
            if centering == 'avg':
                hidens = numpy.where(dens >= dens_limit)[0]
                while hidens.size < nparticles:
                    dens_limit /= 2
                    hidens = numpy.where(dens >= dens_limit)[0]
                if verbose:
                    print ('Center averaged over %d particles' %hidens.size)
                    print ('Center averaged over all particles with density '\
                               'greater than %.2e particles/cc' %dens_limit)
                #Center on highest density clump, rejecting outliers:
                cx = numpy.average(reject_outliers(x[hidens]))
                cy = numpy.average(reject_outliers(y[hidens]))
                cz = numpy.average(reject_outliers(z[hidens]))
                print 'Density averaged box center: %.3e %.3e %.3e' %(cx,cy,cz)
                cu = numpy.average(reject_outliers(u[hidens]))
                cv = numpy.average(reject_outliers(v[hidens]))
                cw = numpy.average(reject_outliers(w[hidens]))
                print 'Density averaged center velocity: %.3e %.3e %.3e' %(cu,cv,cw)
            elif centering == 'max':
                center = dens.argmax()
                cx,cy,cz = x[center], y[center], z[center]
                print 'Density maximum box center: %.3e %.3e %.3e' %(cx,cy,cz)
                cu,cv,cw = u[center], v[center], w[center]
                print 'Density maximum center velocity: %.3e %.3e %.3e' %(cu,cv,cw)
        else:
            raise KeyError("'avg' and 'max' centering require gas density")
    elif centering == 'box':
        cx = (x.max() + x.min())/2
        cy = (y.max() + y.min())/2
        cz = (z.max() + z.min())/2
        print 'Simple box center: %.3e %.3e %.3e' %(cx,cy,cz)
        cu, cv, cw = 0, 0, 0
    else:
        raise KeyError("Centering options are 'avg', 'max', and 'box'")
    return (cx, cy, cz), (cu, cv, cw)

def center_box(x, y, z, u=None, v=None, w=None, center=None, vcenter=None, **kwargs):
    if u is not None or v is not None or w is not None:
        if u is not None and v is not None and w is not None:
            center_velocity = True
        else:
            raise KeyError("Specifying u or v or w requires specifying all three.")
    else:
        center_velocity = False
    dens = kwargs.pop('density', None)
    centering = kwargs.get('centering', None)
    if center:
        cx = center[0]
        cy = center[1]
        cz = center[2]
        if vcenter is not None:
            cu = vcenter[0]
            cv = vcenter[1]
            cw = vcenter[2]
        else:
            cu, cv, cw = 0, 0, 0
    elif centering:
        if center_velocity:
            (cx,cy,cz), (cu,cv,cw) = find_center_vcenter(x,y,z,u,v,w, dens, **kwargs)
        else:
            (cx,cy,cz) = find_center(x,y,z, dens, **kwargs)
    else:
        print "WARNING! NO CENTER OR CENTERING ALGORITHM SPECIFIED!"
        print "Attempting simple box centering..."
        if center_velocity:
            (cx,cy,cz), (cu,cv,cw) = find_center_vcenter(x,y,z,u,v,w, dens, **kwargs)
        else:
            (cx,cy,cz) = find_center(x,y,z, dens, **kwargs)

    x -= cx
    y -= cy
    z -= cz
    if center_velocity:
        u -= cu
        v -= cv
        w -= cw
        return x, y, z, u, v, w
    else:
        return x, y, z

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

# coordinates.py
# Jacob Hummel
import numpy as np

#===============================================================================
def cartesian_to_spherical(x, y, z):
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arccos(z/r)
    phi = np.arctan2(y,x)
    return r,theta,phi

def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def cartesian_to_cylindrical(x ,y, z):
    r = np.sqrt(np.square(x) + np.square(y))
    theta = np.arctan2(y,x)
    return r,theta,z

def cylindrical_to_cartesian(r, theta, z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x,y,z

def define_unit_vectors(pos):
    e1 = np.zeros_like(pos)
    e2 = np.zeros_like(pos)
    e3 = np.zeros_like(pos)
    e1[...,0] = 1.
    e2[...,1] = 1.
    e3[...,2] = 1.
    return e1, e2, e3

def cartesian_to_spherical_unit_vectors(xyz):
    unit_x, unit_y, unit_z = define_unit_vectors(xyz)
    r,theta,phi = cartesian_to_spherical(xyz[...,0],xyz[...,1],xyz[...,2])
    r = r[...,np.newaxis]
    theta = theta[...,np.newaxis]
    phi = phi[..., np.newaxis]
    unit_r = xyz / r
    unit_phi = unit_y * np.cos(phi) - unit_x * np.sin(phi)
    unit_theta = (unit_x * np.cos(theta) * np.cos(phi)
                  + unit_y * np.cos(theta) * np.sin(phi)
                  - unit_z * np.sin(theta))
    return unit_r, unit_theta, unit_phi

def spherical_to_cartesian_unit_vectors(sph):
    unit_r, unit_theta, unit_phi = define_unit_vectors(sph)
    r = sph[...,0, np.newaxis]
    theta = sph[...,1, np.newaxis]
    phi = sph[...,2, np.newaxis]
    unit_x = (unit_r * np.sin(theta) * np.cos(phi)
              + unit_theta * np.cos(theta) * np.cos(phi)
              - unit_phi * np.sin(phi))
    unit_y = (unit_r * np.sin(theta) * np.sin(phi)
              + unit_theta * np.cos(theta) * np.sin(phi)
              + unit_phi * np.cos(phi))
    unit_z = unit_r * np.cos(theta) - unit_theta * np.sin(theta)
    return unit_x, unit_y, unit_z

def cartesian_to_spherical_velocities(xyz, uvw):
    unit_r, unit_theta, unit_phi = cartesian_to_spherical_unit_vectors(xyz)
    # np.einsum('ij,ij->i',...) is einstein notation for dot product.
    vr = np.einsum('ij,ij->i', uvw, unit_r)
    vtheta = np.einsum('ij,ij->i', uvw, unit_theta)
    vphi = np.einsum('ij,ij->i', uvw, unit_phi)
    return vr, vtheta, vphi

def spherical_to_cartesian_velocities(sph, vsph):
    unit_x, unit_y, unit_z = spherical_to_cartesian_unit_vectors(sph)
    # np.einsum('ij,ij->i',...) is einstein notation for dot product.
    u = np.einsum('ij,ij->i', vsph, unit_x)
    v = np.einsum('ij,ij->i', vsph, unit_y)
    w = np.einsum('ij,ij->i', vsph, unit_z)
    return u, v, w

def rotation_matrix(axis, angle):
    if axis == 'x':
        rot = ((1., 0., 0.),
               (0., numpy.cos(angle), -numpy.sin(angle)),
               (0., numpy.sin(angle), numpy.cos(angle)))
    elif axis == 'y':
        rot = ((numpy.cos(angle), 0.,numpy.sin(angle)),
               (0., 1., 0.),
               (-numpy.sin(angle), 0., numpy.cos(angle)))
    elif axis == 'z':
        rot = ((numpy.cos(angle), -numpy.sin(angle), 0.),
               (numpy.sin(angle), numpy.cos(angle), 0.),
               (0., 0.,1.))
    else:
        raise KeyError("{} is not a valid axis choice".format(axis))
    rot = numpy.asarray(rot)
    return rot

def rotate(coords, axis, angle):
    rot = rotation_matrix(axis,angle)
    print "Rotating about the {}-axis by {:6.3f} radians.".format(axis,angle)
    print "Rotation Matrix:"
    print rot
    return numpy.dot(coords,rot)

# coordinates.py
# Jacob Hummel
import numpy as np

#===============================================================================
def cartesian_to_spherical(coords):
    coords['r_sph'] = np.sqrt(np.square(coords.x) + np.square(coords.y) + np.square(coords.z))
    coords['theta_sph'] = np.arccos(coords.z/coords.r_sph)
    coords['phi_sph'] = np.arctan2(coords.y,coords.x)

def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def cartesian_to_cylindrical(coords):
    coords['r_cyl'] = np.sqrt(np.square(coords.x) + np.square(coords.y))
    coords['theta_cyl'] = np.arctan2(coords.y,coords.x)
    #return coords.drop(['x', 'y', 'z'], axis=1)

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

def cartesian_to_spherical_unit_vectors(coords):
    unit_x, unit_y, unit_z = define_unit_vectors(coords)
    r,theta,phi = cartesian_to_spherical(coords[...,0],coords[...,1],coords[...,2])
    r = r[...,np.newaxis]
    theta = theta[...,np.newaxis]
    phi = phi[..., np.newaxis]
    unit_r = coords / r
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

def cartesian_to_spherical_velocities(coords, uvw):
    unit_r, unit_theta, unit_phi = cartesian_to_spherical_unit_vectors(coords)
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
               (0., np.cos(angle), -np.sin(angle)),
               (0., np.sin(angle), np.cos(angle)))
    elif axis == 'y':
        rot = ((np.cos(angle), 0.,np.sin(angle)),
               (0., 1., 0.),
               (-np.sin(angle), 0., np.cos(angle)))
    elif axis == 'z':
        rot = ((np.cos(angle), -np.sin(angle), 0.),
               (np.sin(angle), np.cos(angle), 0.),
               (0., 0.,1.))
    else:
        u = axis[0]
        v = axis[1]
        w = axis[2]
        u2= u*u
        v2= v*v
        w2= w*w
        cos = np.cos(angle)
        sin = np.sin(angle)
        rot = ((u2 + (1-u2)*cos, u*v*(1-cos) - w*sin, u*w*(1-cos) + v*sin),
               (u*v*(1-cos) + w*sin, v2 + (1-v2)*cos, v*w*(1-cos) - u*sin),
               (u*w*(1-cos) - v*sin, v*w*(1-cos) + u*sin, w2 + (1-w2)*cos))
    rot = np.asarray(rot)
    return rot

def rotate(coords, axis, angle, verbose=False):
    rot = rotation_matrix(axis,angle)
    if verbose:
        print "Rotating about the {}-axis by {:6.3f} radians.".format(axis,angle)
        print "Rotation Matrix:"
        print rot
    return np.dot(coords,rot)

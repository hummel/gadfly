# coordinates.py
# Jacob Hummel
import numpy as np

#===============================================================================
def cart2sph(x, y, z):
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arccos(z/r)
    phi = np.arctan2(y,x)
    return r,theta,phi

def sph2cart(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def cart2cyl(x ,y, z):
    r = np.sqrt(np.square(x) + np.square(y))
    theta = np.arctan2(y,x)
    return r,theta,z

def cyl2cart(r, theta, z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x,y,z

def cartesian_unit_vectors(xyz):
    unit_x = np.zeros_like(xyz)
    unit_y = np.zeros_like(xyz)
    unit_z = np.zeros_like(xyz)
    unit_x[...,0] = 1.
    unit_y[...,1] = 1.
    unit_z[...,2] = 1.
    return unit_x, unit_y, unit_z

def spherical_unit_vectors(xyz):
    unit_x, unit_y, unit_z = cartesian_unit_vectors(xyz)
    r,theta,phi = cart2sph(xyz[...,0],xyz[...,1],xyz[...,2])
    r = r[...,np.newaxis]
    theta = theta[...,np.newaxis]
    phi = phi[..., np.newaxis]
    unit_r = xyz / r
    unit_phi = unit_y * np.cos(phi) - unit_x * np.sin(phi)
    unit_theta = (unit_x * np.cos(theta) * np.cos(phi)
                  + unit_y * np.cos(theta) * np.sin(phi)
                  - unit_z * np.sin(theta))
    return unit_r, unit_theta, unit_phi

def sph2cart_unit_vectors(sph):
    r = sph[...,0, np.newaxis]
    theta = sph[...,1, np.newaxis]
    phi = sph[...,2, np.newaxis]
    x,y,z = sph2cart(r, theta, phi)
    xyz = np.column_stack((x,y,z))
    unit_r, unit_theta, unit_phi = spherical_unit_vectors(xyz)
    unit_x = (unit_r * np.sin(theta) * np.cos(phi)
              + unit_theta * np.cos(theta) * np.cos(phi)
              - unit_phi * np.sin(phi))
    unit_y = (unit_r * np.sin(theta) * np.sin(phi)
              + unit_theta * np.cos(theta) * np.sin(phi)
              + unit_phi * np.cos(phi))
    unit_z = unit_r * np.cos(theta) - unit_theta * np.sin(theta)
    return unit_x, unit_y, unit_z

def cart2sph_velocities(xyz, uvw):
    unit_r, unit_theta, unit_phi = spherical_unit_vectors(xyz)
    # np.einsum('ij,ij->i',...) is einstein notation for dot product.
    vr = np.einsum('ij,ij->i', uvw, unit_r)
    vtheta = np.einsum('ij,ij->i', uvw, unit_theta)
    vphi = np.einsum('ij,ij->i', uvw, unit_phi)
    return vr, vtheta, vphi

def sph2cart_velocities(sph, vsph):
    unit_x, unit_y, unit_z = sph2cart_unit_vectors(sph)
    # np.einsum('ij,ij->i',...) is einstein notation for dot product.
    u = np.einsum('ij,ij->i', vsph, unit_x)
    v = np.einsum('ij,ij->i', vsph, unit_y)
    w = np.einsum('ij,ij->i', vsph, unit_z)
    return u, v, w

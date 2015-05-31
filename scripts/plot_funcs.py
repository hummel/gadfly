# plot_funcs.py
# Jacob Hummel
"""
	This file contains convenience functions for common analysis plots.
"""


import plotting
import multiplot

#===============================================================================
def plot_temp(snapshot,wpath=None):
    fig = multiplot.Phase(snapshot)
    fig.plot('temp')
    if wpath:
        fpath = wpath + '/gas/temp/'
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fig.save(fpath+'{:0>4}-temp.png'.format(snapshot.number))

def plot_radial_temp(snapshot,wpath=None):
    fig = multiplot.Phase(snapshot)
    fig.plot('radial-temp')
    if wpath:
        fpath = wpath + '/gas/radial-temp/'
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fig.save(fpath+'{:0>4}-temp.png'.format(snapshot.number))

def plot_gas_fraction(snapshot,wpath=None):
    fig = plotting.Quad(snapshot)
    fig.plot('temp','electron_frac','h2frac','HDfrac')
    if wpath:
        fpath = wpath + '/gas/frac/'
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fig.save(fpath+'{:0>4}-frac.png'.format(snapshot.number))

def disk_density_structure(snapshot, wpath=None):
    fig = plotting.Image(snapshot, track_sinks=True)
    if snapshot.sim.batch_viewscale:
        scale = snapshot.sim.batch_viewscale
    else:
        scale = '5000AU'
    for view in ['xy', 'xz', 'yz']:
        fig.density(scale, view, clim=(8,12), centering='avg')
        if wpath:
            fpath = wpath + '/disk/{}/'.format(view)
            if not os.path.exists(fpath):
                os.makedirs(fpath)
            fig.save(fpath+'{:0>4}-disk-{}.png'.format(snapshot.number, view))

def halo_density_structure(snapshot, wpath=None):
    fig = plotting.Image(snapshot)
    if snapshot.sim.batch_viewscale:
        scale = snapshot.sim.batch_viewscale
    else:
        scale = '1pc'
    fig.density(scale, 'xy', centering='avg')
    if wpath:
        fpath = wpath + '/halo/{}/'.format(scale)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fig.save(fpath+'{:0>4}-halo-{}.png'.format(snapshot.number, scale))

def box_structure(snapshot, wpath=None):
    snapshot.sim.set_coordinate_system('comoving')
    fig = plotting.Image(snapshot)
    if snapshot.sim.batch_viewscale:
        scale = snapshot.sim.batch_viewscale
    else:
        scale = '140kpc'
    fig.density(scale, 'xy', depth=.5, centering='box')
    if wpath:
        fpath = wpath + '/box/'
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        fig.save(fpath+'{:0>4}-box.png'.format(snapshot.number))

def disk_rotation(snapshot, view, rot_axis, n, wpath=None):
    fig = plotting.Image(snapshot, track_sinks=True)
    if snapshot.sim.batch_viewscale:
        scale = snapshot.sim.batch_viewscale
    else:
        scale = '10000AU'
    axis_init = view.get(rot_axis, 0)
    count = 0
    for angle in numpy.linspace(0,2*numpy.pi, n):
        view[rot_axis] = axis_init + angle
        fig.density(scale, view, clim=(7,12), centering='avg')
        if wpath:
            fpath = wpath + '/disk/{}rotation/'.format(rot_axis)
            if not os.path.exists(fpath):
                os.makedirs(fpath)
            fig.save(fpath+'{:0>4}-disk-{}.png'.format(count,rot_axis))
            count += 1


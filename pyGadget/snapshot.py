# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""
import os
import string
import threading
import h5py
import numpy

import hdf5
import nbody
import sph
import plotting

class File:
    """
    Class for Gadget2 HDF5 snapshot files.
    """
    def __init__(self, sim, filename, **kwargs):
        self.sim = sim
        self.filename = filename
        f = filename.replace('.hdf5','')
        self.number = int(f.split('_')[-1])
        self.file_id = h5py.File(filename, 'r')
        self.header = hdf5.Header(self.file_id)
        self.dm = nbody.PartTypeDM(self.file_id, sim.units, **kwargs)
        self.gas = sph.PartTypeSPH(self.file_id, sim.units, **kwargs)
        self.sinks = []
        if kwargs.get('track_sinks', False):
            self.sinks = self.gas.get_sink_properties()

    def update_sink_coordinates(self, x,y,z):
        for s in self.sinks:
            s.update_coordinates(x,y,z)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        self.file_id.close()

#===============================================================================
class Loader(threading.Thread):
    def __init__(self, load_function, file_queue, data_queue):
        self.file_queue = file_queue
        self.data_queue = data_queue
        self.load_function = load_function
        threading.Thread.__init__(self)
        
    def run(self):
        lock = threading.Lock()
        while 1:
            args = self.file_queue.get()
            if args is None:
                self.data_queue.put(None)
                break # reached end of queue
            fname = args[0]
            lock.acquire()
            print 'loading snapshot', fname
            lock.release()
            try:
                snapshot = self.load_function(*args)
                self.data_queue.put(snapshot)
            except IOError:
                lock.acquire()
                print 'Warning: snapshot '+str(fname)+' not found!'
                lock.release()
                pass

#===============================================================================
def plot_temp(snapshot,wpath=None):
    fig = plotting.Phase(snapshot)
    fig.plot('temp')
    if wpath:
        fpath = wpath + '/gas/temp/'
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
        fig.density(scale, view, clim=(7,12), centering='avg')
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

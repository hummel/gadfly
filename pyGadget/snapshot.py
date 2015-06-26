# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""
import os
import string
import threading
import Queue
import h5py
import numpy

import hdf5
import nbody
import sph

class File(object):
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

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['file_id']
        return result

#    def __setstate__(self, in_dict):
#        self.__dict__ = in_dict
#        self.file_id = h5py.File(self.filename, 'r')
#        self.header = hdf5.Header(self.file_id)

    def update_sink_coordinates(self, x,y,z):
        for s in self.sinks:
            s.update_coordinates(x,y,z)

    def update_sink_velocities(self, u,v,w):
        for s in self.sinks:
            s.update_velocities(u,v,w)

    def update_sink_frame_ofR(self, xyz, uvw):
        self.update_sink_coordinates(xyz[:,0], xyz[:,1], xyz[:,2])
        self.update_sink_velocities(uvw[:,0], uvw[:,1], uvw[:,2])
        
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
            try:
                args = self.file_queue.get(timeout=1)
            except Queue.Empty:
                break # reached end of queue
            if args is None:
                self.data_queue.put(None)
            else:
                fname = args[0]
                lock.acquire()
                print 'loading snapshot', fname
                lock.release()
                try:
                    snapshot = self.load_function(*args)
                    snapshot.close()
                    self.data_queue.put(snapshot)
                except IOError:
                    lock.acquire()
                    print 'Warning: snapshot '+str(fname)+' not found!'
                    lock.release()
                    pass


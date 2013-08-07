# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""
import string
import threading
import h5py

import hdf5
import nbody
import sph
import plotting

class File:
    """
    Class for Gadget2 HDF5 snapshot files.
    """
    def __init__(self, filename):
        self.filename = filename
        f = filename.replace('.hdf5','')
        self.number = int(f.split('_')[-1])
        self.file_id = h5py.File(filename, 'r')
        self.header = hdf5.Header(self.file_id)
        self.dm = nbody.PartTypeDM(self.file_id)
        self.gas = sph.PartTypeSPH(self.file_id)
        
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
            print 'loading ' + fname
            lock.release()
            try:
                snapshot = self.load_function(*args)
                self.data_queue.put(snapshot)
            except IOError:
                lock.acquire()
                print 'Warning: '+fname+' not found!'
                lock.release()
                pass

#===============================================================================
def plot_temp(snapshot,wpath=None):
    fig = plotting.Phase(snapshot)
    fig.plot('temp')
    if wpath:
        fig.save(wpath+'/gas/temp/{:0>4}-temp.png'.format(snapshot.number))

def plot_gas_fraction(snapshot,wpath=None):
    fig = plotting.Quad(snapshot)
    fig.plot('temp','electron_frac','h2frac','HDfrac')
    if wpath:
        fig.save(wpath+'/gas/frac/{:0>4}-frac.png'.format(snapshot.number))

# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""
import threading
import h5py
import hdf5
import nbody
import sph

class File:
    """
    Class for Gadget2 HDF5 snapshot files.
    """
    def __init__(self, filename):
        self.filename = filename
        self.number = int(filename[-8:-5])
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
        while 1:
            args = self.file_queue.get()
            if args is None:
                self.data_queue.put(None)
                break # reached end of queue
            fname = args[0]
            print 'loading ' + fname
            try:
                snapshot = self.load_function(*args)
                self.data_queue.put(snapshot)
            except IOError:
                print 'Warning: '+fname+' not found!'
                pass


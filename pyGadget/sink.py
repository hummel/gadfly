# sinks.py
# Jacob Hummel
"""
Classes and routines for analyzing sink data output by gadget.
"""
import sys
import warnings
import numpy as np
from astropy.io import ascii
import units
#===============================================================================

class Sink(object):
    def __init__(self,**properties):
        super(Sink,self).__init__()
        self.mass = properties.pop('m', None)
        self.x = properties.pop('x', None)
        self.y = properties.pop('y', None)
        self.z = properties.pop('z', None)
        self.pos = properties.pop('pos', None)
        self.radius = properties.pop('r', None)
        self.energy = properties.pop('e', None)
        self.pressure = properties.pop('p', None)
        self.npart_acc = properties.pop('n', None)
        self.pid = properties.pop('pid', None)
        self.index = properties.pop('index', None)

    def update_coordinates(self, x,y,z):
        self.x = x[self.index]
        self.y = y[self.index]
        self.z = z[self.index]

class SinkData(object):
    def __init__(self,path):
        super(SinkData,self).__init__()
        ### Read in the data
        try:
            sinkdata = ascii.read(path+'/sinkdat')
        except IOError:
            raise IOError("Specified sinkmasses file not found!")
        try:
            sinkmasses = ascii.read(path+'/sinkmasses')
        except IOError:
            raise IOError("Specified sinkmasses file not found!")

        self.time = sinkdata['col1']
        self.npart_acc = sinkdata['col2']
        self.radius = sinkdata['col3']
        self.part_internal_energy = sinkdata['col4']
        self.entropy = sinkdata['col5']
        self.part_id = sinkdata['col6']
        self.sink_id  = sinkdata['col7']
        self.pressure = sinkdata['col8']
        self.a = self.time # Scale Facor

        # Restrict to real sinks
        IDs = np.unique(sinkmasses['col2'])
        real = np.in1d(self.sink_id, IDs)
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][real]

        h = 0.7 #Hubble Parameter
        h2 = h*h
        a3 = self.a**3
        ### Convert units
        self.time = self.time*units.Time_yr
        # npart_acc is a simple integer (no units)
        self.pressure = self.pressure*units.Pressure_cgs*h2/(a3**1.4)
        self.radius = self.radius*units.Length_AU*h

        good = np.where(self.radius > 10)[0]
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][good]

class SingleSink(SinkData):
    '''
    Select sink data for a single sink particle.

    possible options:
    nform: select the n(th) sink to form.
    ID: select sink by ID.
    '''
    def __init__(self, path, nform=None, id_=None):
        super(SingleSink,self).__init__(path)
        unique = np.unique(self.sink_id)

        if((nform is None) and (id_ is None)):
            print "No sink specified: Selecting first sink to form..."
            nform = 1
        if nform:
            print "Key set: nform =", nform
            if id_: warnings.warn("nform key overrides id_")
            # Select n(th) sink to form
            new = []
            i = 0
            while len(new) < nform:
                if self.sink_id[i] not in new:
                    new.append(self.sink_id[i])
                i += 1
            id_ = new[-1]
        elif id_ is not None:
            print "Key set: id_ =", id_
        else:
            raise RuntimeError("Execution should not have reached this point!")
        print "Using sink ID", id_
        
        # Restrict to a single sink
        lines = np.where(self.sink_id == id_)
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][lines]

        # Select final output for each timestep
        tsteps = np.unique(self.time)
        selection = []
        for t in tsteps:
            times = np.where(self.time == t)[0]
            selection.append(times[-1])
        for key in vars(self).keys():
                vars(self)[key] = vars(self)[key][selection]
        self.sink_id = id_

        # Calculate sink mass at each timestep
        self.mass = np.zeros_like(self.time)
        for i in xrange(self.mass.size):
            self.mass[i] = 0.015*self.npart_acc[:i+1].sum()

        # Finally, record total number of sinks found.
        self.all_ids = unique

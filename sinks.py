# sinks.py
# Jacob Hummel
"""
Classes and routines for analyzing sink data output by gadget.
"""
import numpy as np
import asciitable
#===============================================================================

class SinkData(object):
    def __init__(self,filename):
        super(SinkData,self).__init__()
        ### Read in the data
        try:
            sinkdata = asciitable.read(filename)
        except IOError:
            print "Specified sinkdata file not found!"

        self.time = sinkdata['col1']
        self.npart_acc = sinkdata['col2']
        self.r_sink_phys = sinkdata['col3']
        self.part_internal_energy = sinkdata['col4']
        self.sink_entropy = sinkdata['col5']
        self.part_id = sinkdata['col6']
        self.sink_id  = sinkdata['col7']
        self.sink_pressure = sinkdata['col8']

class SingleSink(SinkData):
    def __init__(self, filename, nform=None, ID=None):
        super(SingleSink,self).__init__(filename)
        if id:
            # Restrict to a single sink
            lines = np.where(self.sink_id == ID)
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
        '''
        if number <= self.unique_id.size:
            self.sink_num = number
        else:
            raise KeyError
        '''

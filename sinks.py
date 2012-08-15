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
    def __init__(self, filename, number, sink_id=False):
        super(SingleSink,self).__init__(filename)
        self.unique_id, unique_index = np.unique(self.sink_id,
                                                      return_index=True)
        if number <= self.unique_id.size:
            self.sink_num = number
        else:
            raise KeyError
        


        '''

    # Select only main sink
    sink0 = sinkdata['col7'][0]
    sinkdata = sinkdata[np.where(sinkdata['col7'] == sink0)]
    # Select final output for each timestep
    tpoints = sinkdata['col1']
    unique = np.unique(tpoints)
    selection = []
    for t in unique:
        times = np.where(tpoints == t)[0]
        selection.append(times[-1])
    sinkdata = sinkdata[selection]
    print sinkdata
        '''

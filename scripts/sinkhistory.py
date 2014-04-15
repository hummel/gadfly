#!/usr/bin/env python
# sinkhistory.py
# Jacob Hummel

import sys
import pandas as pd
import pyGadget
#===============================================================================
def get_sinkdata(key, sinkdata):
    snap = sim.load_snapshot(key)
    print "Snapshot {}: z={:.3f}".format(key, snap.header.Redshift)
    t = snap.header.Time
    for s in snap.sinks:
        try:
            sinkdata[s.pid].append((t, s.mass, s.x, s.y, s.z, s.vx, s.vy, s.vz))
        except KeyError:
            sinkdata[s.pid] = [(t, s.mass, s.x, s.y, s.z, s.vx, s.vy, s.vz)]
    snap.gas.cleanup()
    snap.close()
    return sinkdata

if __name__ == '__main__':
    simname = sys.argv[1]
    sim = pyGadget.sim.Simulation(simname, track_sinks=True)
    keys = sim.snapfiles.keys()
    keys.sort()
    sinkdata = {}
    for key in keys:
        sinkdata = get_sinkdata(key, sinkdata)

    data = {}
    for key in sinkdata.keys():
        data[key] = pd.DataFrame(sinkdata[key], columns=('time', 'mass',
                                                         'x', 'y', 'z',
                                                         'u', 'v', 'w'))
    sd = pd.Panel(data)
    store = pd.HDFStore(sim.plotpath+'sinkdata.hdf5')
    store[sim.name] = sd
    store.close()

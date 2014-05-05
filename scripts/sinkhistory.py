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
    t = snap.header.Time * pyGadget.units.Time_yr - sim.tsink
    for s in snap.sinks:
        sinkdata.append((s.pid, t, s.mass, s.x, s.y, s.z, s.vx, s.vy, s.vz))
    snap.gas.cleanup()
    snap.close()
    return sinkdata

if __name__ == '__main__':
    simname = sys.argv[1]
    sim = pyGadget.sim.Simulation(simname, track_sinks=True)
    keys = sim.snapfiles.keys()
    keys.sort()
    sinkdata = []
    for key in keys:
        sinkdata = get_sinkdata(key, sinkdata)

    data = pd.DataFrame(sinkdata, columns=('ID', 'time', 'mass', 'x', 'y', 'z',
                                           'u', 'v', 'w'))

    store = pd.HDFStore(sim.plotpath+'sinkdata.hdf5')
    store[sim.name] = data
    store.close()

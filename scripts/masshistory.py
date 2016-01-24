#!/usr/bin/env python
# sinkhistory.py
# Jacob Hummel

import sys
import numpy as np
import pandas as pd
import pyGadget
#===============================================================================
def get_cdgm(sim, tracer_pid, key):
    snap = sim.load_snapshot(key, 'particleIDs', 'masses', 'coordinates', 'velocities',
                             'ndensity')
    z = snap.header.Redshift
    if sim.tsink:
        t = snap.header.Time * pyGadget.units.Time_yr - sim.tsink
    else:
        t = snap.header.Time * pyGadget.units.Time_yr
    if t > 0:
        print "Snapshot {}: z={:.3f}".format(key, z),
        print "t_sink={:.1f}".format(t)
    else:
        print "Snapshot {}: z={:.3f}".format(key, z)
    pids = snap.gas.get_PIDs()
    tracer_idx = np.where(pids == tracer_pid)[0]
    del pids

    xyz = snap.gas.get_coords(unit='pc', system='cartesian')
    center = xyz[tracer_idx][0]
    tracer_pos = (center[0], center[1], center[2])
    pos = snap.gas.get_coords(unit='pc', system='spherical', center=tracer_pos)
    r = pos[:,0]
    del xyz, pos

    dens = snap.gas.get_number_density()
    mass = snap.gas.get_masses()
    mdata = [z,t]
    for dlim in [10,100,1e4,1e6,1e8,1e9,1e10,1e11]:
        (m,) = pyGadget.analyze.data_slice(dens >= dlim, mass)
        mdata.append(m.sum())
    for rlim in [1000, 100, 10, 1, .1, 0.04848, 0.02424, 0.004848]:
        (m,) = pyGadget.analyze.data_slice(r <= rlim, mass)
        mdata.append(m.sum())
    snap.gas.cleanup()
    snap.close()
    return mdata

if __name__ == '__main__':
    simname = sys.argv[1]
    tracer_id = int(sys.argv[2])
    sim = pyGadget.sim.Simulation(simname, length='pc', track_sinks=True)
    keys = sim.snapfiles.keys()
    keys.sort()
    ikeys = keys
    nkeys = ['10cc', '100cc', '1e4cc', '1e6cc', '1e8cc', '1e9cc', '1e10cc', '1e11cc']
    rkeys = ['1kpc', '100pc', '10pc', '1pc', '.1pc', '1e4AU', '5e3AU', '1e3AU']
    col_names = ['z', 'time']
    col_names.extend(nkeys)
    col_names.extend(rkeys)
    df = pd.DataFrame(index=ikeys, columns=col_names)
    for key in ikeys:
        try:
            mdata = get_cdgm(sim, tracer_id, key)
            df.loc[key] = mdata
        except(IOError):
            pass

    store = pd.HDFStore(sim.plotpath+'mass_history.hdf5')
    store[sim.name] = df
    store.close()

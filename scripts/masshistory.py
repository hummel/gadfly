#!/usr/bin/env python
# sinkhistory.py
# Jacob Hummel

import sys
import pandas as pd
import pyGadget
#===============================================================================
def get_cdgm(sim, key):
    snap = sim.load_snapshot(key, 'masses', 'coordinates', 'velocities',
                             'ndensity')
    z = snap.header.Redshift
    t = snap.header.Time * pyGadget.units.Time_yr - sim.tsink
    if t > 0:
        print "Snapshot {}: z={:.3f}".format(key, z),
        print "t_sink={:.1f}".format(t)
    else:
        print "Snapshot {}: z={:.3f}".format(key, z)
    pos = snap.gas.get_coords(unit='pc', system='spherical', centering='avg')
    r = pos[:,0]
    del pos
    dens = snap.gas.get_number_density()
    mass = snap.gas.get_masses()
    mdata = [z,t]
    for dlim in [10,100,1e4,1e8,1e10]:
        (m,) = pyGadget.analyze.data_slice(dens >= dlim, mass)
        mdata.append(m.sum())
    for rlim in [100, 10, 1, .1, 0.04848, 0.02424, 0.004848]:
        (m,) = pyGadget.analyze.data_slice(r <= rlim, mass)
        mdata.append(m.sum())
    snap.gas.cleanup()
    snap.close()
    return mdata

if __name__ == '__main__':
    simname = sys.argv[1]
    sim = pyGadget.sim.Simulation(simname, length='pc', track_sinks=True)
    keys = sim.snapfiles.keys()
    keys.sort()
#    ikeys = keys[:500] + keys[500::2]
    ikeys = keys
    nkeys = ['10cc', '100cc', '1e4cc', '1e8cc', '1e10cc']
    rkeys = ['100pc', '10pc', '1pc', '.1pc', '1e4AU', '5e3AU', '1e3AU']
    col_names = ['z', 'time']
    col_names.extend(nkeys)
    col_names.extend(rkeys)
    df = pd.DataFrame(index=ikeys, columns=col_names)
    for key in ikeys:
        try:
            mdata = get_cdgm(sim, key)
            df.loc[key] = mdata
        except(IOError):
            pass

    store = pd.HDFStore(sim.plotpath+'mass_history.hdf5')
    store[sim.name] = df
    store.close()

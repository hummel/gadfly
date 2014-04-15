import sys
import pyGadget

if __name__ == '__main__':
    simname = sys.argv[1]
    sim = pyGadget.sim.Simulation(simname, refine=False, track_sinks=True)
    keys = sim.snapfiles.keys()
    keys.sort()
    sim.track_sinks(False)
    for key in keys:
        snap = sim.load_snapshot(key)
        print "Snapshot {}: z={:.3f}".format(key, snap.header.Redshift)
        if snap.header.ScaleFactor > sim.sink1.a[0]:
            for sink in sim.sinks:
                if (snap.header.ScaleFactor >= sink.a[0] and
                    snap.header.ScaleFactor <= sink.a[-1]):
                    disk = pyGadget.sink.AccretionDisk(sim, sink)
                    try:
                        disk.load(snap.number, verbose=True)
                    except IndexError:
                        pass
        snap.gas.cleanup()
        snap.close()

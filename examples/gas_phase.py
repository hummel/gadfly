#!/usr/bin/env python
# gas_phase.py
# Jacob Hummel

import sys
import pyGadget

#===============================================================================
if __name__ == '__main__':
    if ((len(sys.argv) not in [3,4,5]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gas_phase.py [key] [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        print '   Option 2: python gas_phase.py [key] [simulation name] '\
            '[single snapshot]'
        print '   Option 3: python gas_phase.py [key] [simulation name] '\
            '(creates a plot for every snapshot)'
        sys.exit()

    key = sys.argv[1]
    if key == 'temp':
        plot_func = pyGadget.snapshot.plot_temp
        data = ['ndensity','temp']
    elif key == 'radial-temp':
        plot_func = pyGadget.snapshot.plot_radial_temp
        data = ['ndensity','temp','coordinates']
    elif key == 'frac':
        plot_func = pyGadget.snapshot.plot_gas_fraction
        data = ['ndensity','temp','electron_frac','h2frac','HDfrac']
    else:
        raise KeyError
    simname = sys.argv[2]
    sim = pyGadget.sim.Simulation(simname, length='pc')
                 
    if len(sys.argv) == 4:
        snap = int(sys.argv[3])
        snapshot = sim.load_snapshot(snap,*data)
        plot_func(snapshot, sim.plotpath+sim.name)

    elif len(sys.argv) == 5:
        start = int(sys.argv[3])
        stop = int(sys.argv[4])
        snaps = range(start,stop+1)
        sim.set_snapshots(*snaps)
        sim.multitask(plot_func,*data)

    else:
        sim.multitask(plot_func,*data)



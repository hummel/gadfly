#!/usr/bin/env python
# structure.py
# Jacob Hummel

import sys
import pyGadget

#===============================================================================
help_message = \
    """
Usage::
    Option 1: python gas_phase.py [key] [scale] [simulation name]
              [beginning snapshot] [final snapshot]              
    Option 2: python gas_phase.py [key] [scale] [simulation name]
              [single snapshot]
    Option 3: python gas_phase.py [key] [scale] [simulation name] 
              (creates a plot for every snapshot)

Key Options: 'disk' 'halo' 'box'
Scale: Must be specified as alphanumeric string (no scientific notation)
       unit options: 'AU' 'pc' 'kpc' 'cm'
    """
if __name__ == '__main__':
    if ((len(sys.argv) not in [4,5,6]) or (sys.argv[1] == '-h')):
        print help_message
        sys.exit()
    
    key = sys.argv[1]
    if key == 'disk':
        plot_func = pyGadget.snapshot.disk_density_structure
    elif 'halo' in key:
        plot_func = pyGadget.snapshot.halo_density_structure
    elif key == 'box':
        plot_func = pyGadget.snapshot.box_structure
    else:
        raise KeyError

    scale = sys.argv[2]
    unit = "".join(ch if not ch.isdigit() else "" for ch in scale)
    simname = sys.argv[3]
    sim = pyGadget.sim.Simulation(simname, length=unit)
    sim.set_batch_viewscale(scale)
    if key == 'disk':
        sim.track_sinks()
    if key == 'box':
        sim.refine_by_mass(False)
        sim.set_coordinate_system('comoving')
                 
    data = ['ndensity','coordinates', 'smoothing_length']
    if len(sys.argv) == 5:
        snap = int(sys.argv[4])
        snapshot = sim.load_snapshot(snap,*data)
        plot_func(snapshot, sim.plotpath+sim.name)

    elif len(sys.argv) == 6:
        start = int(sys.argv[4])
        stop = int(sys.argv[5])
        snaps = range(start,stop+1)
        sim.set_snapshots(*snaps)
        sim.multitask(plot_func, *data, parallel=False)

    else:
        sim.multitask(plot_func, *data, parallel=False)



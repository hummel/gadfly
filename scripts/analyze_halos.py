# halos.py
# Jacob Hummel
# quick script for updating all halo databases.
import pyGadget

for simulation in ['stampede/vanilla', 'stampede/XR_sfr_1e-1',
                   'stampede/XR_sfr_1e-2', 'stampede/XR_sfr_1e-3']:
    sim = pyGadget.sim.Simulation(simulation)
    halo = pyGadget.halo.Halo(sim)
    for snap in sim.snapfiles.keys():
        halo.populate(snap)

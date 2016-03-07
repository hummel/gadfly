# halos.py
# Jacob Hummel
# quick script for updating all halo databases.
import sys
import pyGadget

simulation = sys.argv[1]
sim = pyGadget.sim.Simulation(simulation)
halo = pyGadget.halo.Halo(sim)
keys = sim.snapfiles.keys()
keys.sort()
for snap in keys:
    halo.load(snap)

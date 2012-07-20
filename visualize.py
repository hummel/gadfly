# particles3D.py
# plot particle positions in 3D
# Jacob Hummel

import os

import numpy as np
import h5py
from mayavi import mlab

from gadgetHDF5 import *
#===============================================================================
mlab.figure(1, bgcolor=(0, 0, 0), size=(1000, 1000))
mlab.clf()

home = os.getenv('HOME')
datafile = home+'/sim/snapshot_037.hdf5'

data = h5py.File(datafile,'r')
pos = data['PartType1'].get('Coordinates')[...]
mass = data['PartType1'].get('Masses')[...]
density = data['PartType0'].get('Density')[...]
data.close()

stride = 10000
pos = pos[0:-1:stride]
mass = mass[0:-1:stride]*1e11
density = density[0:-1:stride]*1e9


x = pos[:,0]
y = pos[:,1]
z = pos[:,2]

mmin = mass.min()
sizes = np.where(mass == mmin, 1, mass)
sizes = np.where(sizes == 8.*mmin, 2, sizes)
sizes = np.where(sizes == 64.*mmin, 3, sizes)
sizes = np.where(sizes == 512.*mmin, 4, sizes)

#pts = mlab.points3d(x, y, z, sizes, scale_factor=.25, resolution=12)
#mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts))

#source = mlab.pipeline.scalar_scatter(x,y,z)
'''
min = density.min()
max = density.max()
vol = mlab.pipeline.volume(source,
                           vmin=min+0.65*(max-min),
                           vmax=min+0.9*(max-min))
'''
#mlab.view(132, 54, 45, [21, 20, 21.5])

#mlab.show()


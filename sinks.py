# particles3D.py
# plot particle positions in 3D
# Jacob Hummel

import os
import sys
from matplotlib import pyplot

import pyGadget
#===============================================================================

if len(sys.argv) < 2:
    print 'Usage: python sinks.py (simulation name)'
    sys.exit()

simulation = sys.argv[1]
path = os.getenv('HOME')+'/data/sinks/'+simulation+'/sinkdat'
wdir = os.getenv('HOME')+'/data/sinks/'+simulation+'/'
sink = pyGadget.sinks.SingleSink(path,1)
sink.time = sink.time - sink.time[0] # set t=0 at sink formation.

### Create Plot!
print 'Plotting...'
fig = pyplot.figure(1,(12,10))
fig.clf()
ax0 = fig.add_subplot(221)
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(224)

# Particles accreted
ax0.plot(sink.time, sink.npart_acc,'o')
ax0.set_xlabel('Time Since Formation [yr]')
ax0.set_ylabel('Particles Accreted')

# Sink Physical Radius
ax1.plot(sink.time, sink.radius)
ax1.set_xlabel('Time Since Formation [yr]')
ax1.set_ylabel('Sink Radius [AU/h]')

# Sink Entropy
ax2.plot(sink.time, sink.entropy)
ax2.set_xlabel('Time Since Formation [yr]')
ax2.set_ylabel('Sink Entropy [???]')

# Sink Pressure
ax3.plot(sink.time, sink.pressure)
ax3.set_xlabel('Time Since Formation [yr]')
ax3.set_ylabel('Sink Pressure [dyne cm$^{-2}$]')
#ax3.set_yscale('log')
ax3.set_ylim(.1,.5)

fig.subplots_adjust(top=0.94, left=0.1, right=.915)
pyplot.savefig(wdir+'sinkdat.png', 
               bbox_inches='tight')

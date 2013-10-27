#!/usr/bin/env python
# sinks.py
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
path = os.getenv('HOME')+'/data/sinks/'+simulation+'/'
wdir = os.getenv('HOME')+'/data/sinks/'+simulation+'/'
sink1 = pyGadget.sink.SingleSink(path,1)
sink2 = pyGadget.sink.SingleSink(path,2)
t0 = sink1.time[0] # set t=0 at first sink formation.
sink1.time = sink1.time - t0
sink2.time = sink2.time - t0

### Create Plot!
print 'Plotting...'
fig = pyplot.figure(1,(12,10))
fig.clf()
ax0 = fig.add_subplot(221)
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(224)

# Particles accreted
ax0.plot(sink1.time, sink1.npart_acc,'o',label='Sink 1')
ax0.plot(sink2.time, sink2.npart_acc,'o',label='Sink 2')
ax0.set_xlabel('Time Since Formation [yr]')
ax0.set_ylabel('Particles Accreted')
ax0.legend()

# Sink Physical Radius
ax1.plot(sink1.time, sink1.radius)
ax1.plot(sink2.time, sink2.radius)
ax1.set_xlabel('Time Since Formation [yr]')
ax1.set_ylabel('Sink Radius [AU]')

# Sink Mass
ax2.plot(sink1.time, sink1.mass)
ax2.plot(sink2.time, sink2.mass)
ax2.set_xlabel('Time Since Formation [yr]')
ax2.set_ylabel('Sink Mass [M$_{\odot}$]')

# Sink Pressure
ax3.plot(sink1.time, sink1.pressure)
ax3.plot(sink2.time, sink2.pressure)
ax3.set_xlabel('Time Since Formation [yr]')
ax3.set_ylabel('Sink Pressure [dyne cm$^{-2}$]')
#ax3.set_yscale('log')
ax3.set_ylim(.1,.5)

z = 1/sink1.a[0] -1
fig.suptitle('First Sink Formation: z = %.2f' %z)
fig.subplots_adjust(top=0.94, left=0.1, right=.915)
pyplot.savefig(wdir+'sinks.png', 
               bbox_inches='tight')


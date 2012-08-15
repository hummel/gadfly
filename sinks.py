# sinks.py
# Jacob Hummel
"""
Classes and routines for analyzing sink data output by gadget.
"""
import os
import numpy as np
import asciitable
from matplotlib import pyplot
#===============================================================================

def plot_sinkdat(path, write_dir):
    
    # Code unit conversions
    unitMass_g = 1.989e43 # g
    unitLength_cm = 3.085678e21 # cm
    unitVelocity_cgs= 1.0e5
    unitDensity_cgs= unitMass_g / unitLength_cm**3
    unitTime_s= unitLength_cm / unitVelocity_cgs
    unitDensity_cgs= unitMass_g / unitLength_cm**3
    unitPressure_cgs= unitMass_g / unitLength_cm/ unitTime_s**2
    unitEnergy_cgs= unitMass_g * unitLength_cm**2 / unitTime_s**2
    # Fundamental constants
    k_B = 1.3806e-16 # erg/K
    m_H = 1.6726e-24 # g
    GRAVITY = 6.6726e-8 # dyne * cm**2 / g**2
    G = GRAVITY / unitLength_cm**3 * unitMass_g * unitTime_s**2
    X_h = 0.76 # Hydrogen Mass Fraction

    ### Read in the data
    sinkdata = asciitable.read(path)
    # Select only main sink
    sink0 = sinkdata['col7'][0]
    sinkdata = sinkdata[np.where(sinkdata['col7'] == sink0)]
    # Select final output for each timestep
    tpoints = sinkdata['col1']
    unique = np.unique(tpoints)
    selection = []
    for t in unique:
        times = np.where(tpoints == t)[0]
        selection.append(times[-1])
    sinkdata = sinkdata[selection]
    print sinkdata
        
    time = sinkdata['col1']
    npart_acc = sinkdata['col2']
    r_sink_phys = sinkdata['col3']
    part_internal_energy = sinkdata['col4']
    sink_entropy = sinkdata['col5']
    part_id = sinkdata['col6']
    sink_id  = sinkdata['col7']
    sink_pressure = sinkdata['col8']

    print 'Time Range:', time[0], time[-1]

    ### Create Plot!
    ### All plots vs density
    print 'Plotting...'
    fig = pyplot.figure(1,(12,10))
    fig.clf()
    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)
    
    xlims = (0.03825409, 0.03825415)

    # Particles accreted
    ax0.plot(time, npart_acc)
    #ax0.set_xscale('log')
    #ax0.set_yscale('log')

    #ax0.set_xlim(xlims)
    #ax0.set_ylim(7, 1e4)
    ax0.set_xlabel('Time')
    ax0.set_ylabel('Particles Accreted')

    # Sink Physical Radius
    ax1.plot(time, r_sink_phys)
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')

    #ax1.set_xlim(xlims)
    #ax1.set_ylim(1e-12, 1e-2)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Sink Radius [kpc]')

    # Sink Entropy
    ax2.plot(time, sink_entropy)
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')

    #ax2.set_xlim(xlims)
    #ax2.set_ylim(1e-7,1)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Sink Entropy')

    # Sink Pressure
    ax3.plot(time, sink_pressure)
    ax3.set_xscale('log')

    #ax3.set_xlim(xlims)
    #ax3.set_ylim(1,2)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Sink Pressure')

    #title = fig.suptitle('Redshift: %.3f' %(redshift,))
    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    pyplot.savefig(write_dir+'sinkdat.png', 
                   #bbox_extra_artists=(title,),
                   bbox_inches='tight')
#===============================================================================

if __name__ == '__main__':
    #pyplot.ioff()
    wdir = os.getenv('HOME')+'/data/simplots/vanilla-100/'
    path = os.getenv('HOME')+'/research/work/v100/sinkdat'
    plot_sinkdat(path, wdir)

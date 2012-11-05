# gas_temp.py
# Jacob Hummel

import os
import sys
import glob
import numpy
import Queue
import threading
import multiprocessing
from matplotlib import pyplot
import pyGadget

#===============================================================================
def load_snapshot(path):
    snapshot = pyGadget.snapshot.load(path)
    snapshot.gas.load_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.load_internal_energy()
    snapshot.gas.load_gamma()
    snapshot.gas.load_H2_fraction()
    snapshot.close()
    return snapshot

def plot_temp(snapshot,wpath):
    h = snapshot.header.HubbleParam
    redshift = snapshot.header.Redshift
    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()

    # Refine
    minimum = numpy.amin(particle_mass)
    stride = 100
    refined = numpy.where(particle_mass <= minimum)[0]#[0:-1:stride]
    dens = dens[refined]
    temp = temp[refined]

    # Plot!
    fig = pyplot.figure(figsize=(15,10))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.scatter(dens, temp, s=.75, c='k', linewidths=0.0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e13)
    ax.set_ylim(10, 2e4)

    pyplot.title('Redshift: %.2f' %(redshift,))
    ax.set_xlabel('n [cm^-3]')
    ax.set_ylabel('Temperature [K]')
    pyplot.savefig(wpath+'-temp.png',
                   bbox_inches='tight')

#===============================================================================
class Loader(threading.Thread):
    def __init__(self, file_queue, data_queue):
        self.file_queue = file_queue
        self.data_queue = data_queue
        threading.Thread.__init__(self)
        
    def run(self):
        while 1:
            fname = self.file_queue.get()
            if fname is None:
                self.data_queue.put(None)
                break # reached end of queue
            print 'loading', fname
            snapshot = load_snapshot(fname)
            self.data_queue.put(snapshot)

class Worker(threading.Thread):
    def __init__(self, queue,wdir):
        self.__queue = queue
        self.wdir = wdir
        self.procs = []
        threading.Thread.__init__(self)
        
    def run(self):
        while 1:
            snapshot = self.__queue.get()
            if snapshot is None:
                for p in self.procs:
                    p.join()
                break # reached end of queue
            print 'Plotting', snapshot.filename
            procs = []
            p = multiprocessing.Process(target=plot_temp, 
                                        args=(snapshot,
                                              self.wdir+str(snapshot.number)))
            self.procs.append(p)
            p.start()
            
#===============================================================================
def multitask(path,write_dir,start,stop):
    file_queue = Queue.Queue(2)
    data_queue = Queue.Queue(8)
    Loader(file_queue,data_queue).start()
    Worker(data_queue,write_dir).start()
    for snap in xrange(start,stop):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        file_queue.put(fname)
    file_queue.put(None)

#===============================================================================
if __name__ == '__main__':
    pyplot.ioff()
    if ((len(sys.argv) not in [2,3,4]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gas_temp.py [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        print '   Option 2: python gas_temp.py [simulation name] '\
            '[single snapshot]'
        print '   Option 3: python gas_temp.py [simulation name] '\
            '(creates a plot for every snapshot)'
        sys.exit()

    simulation = sys.argv[1]
    path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
    write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'

    if len(sys.argv) == 3:
        snap = sys.argv[2]
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        snapshot = load_snapshot(fname)
        plot_temp(snapshot,write_dir+snap)
        
    elif len(sys.argv) == 4:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])+1
        multitask(path,write_dir,start,stop)

    else:
        files = glob.glob(path+'???.hdf5')
        files.sort()
        start = int(files[0][-8:-5])
        stop = start + len(files)
        multitask(path,write_dir,start,stop)


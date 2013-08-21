# sim.py
# Jacob Hummel
import os
import glob
import Queue
import subprocess
import multiprocessing as mp

import units
import sinks
import snapshot
import plotting

class Simulation(object):
    """
    Class for simulations.  Primarily for gathering metadata.
    """
    def __init__(self, sim_name, **kwargs):
        super(Simulation,self).__init__()
        self.name = sim_name
        self.datapath = kwargs.pop('datapath',os.getenv('HOME')+'/sim/')
        self.plotpath = kwargs.pop('plotpath',(os.getenv('HOME')
                                               +'/data/simplots/'))
        self.sinkpath = kwargs.pop('sinkpath',os.getenv('HOME')+'/data/sinks/')
        
        self.length_unit = kwargs.pop('length_unit',units.Length_AU)
        self.mass_unit = kwargs.pop('mass_unit',units.Mass_sun)
        self.time_unit = kwargs.pop('length_unit',units.Time_yr)
        self.units = {'length':self.length_unit,
                      'mass':self.mass_unit,
                      'time':self.time_unit}
        self.snapfiles = self.find_snapshots()
        try:
            self.sink0 = sinks.SingleSink(self.sinkpath)
            print "Found sinkfiles.  Loading sinkdata."
            self.tform = self.sink0.time[0]
        except IOError:
            self.tform = None

    def find_snapshots(self, *nums):
        path = self.datapath + self.name
        files0 = glob.glob(path+'/snapshot_???.hdf5')
        files1 = glob.glob(path+'/snapshot_????.hdf5')
        files0.sort()
        files1.sort()
        files = files0 + files1
        snapfiles = {}
        for f in files:
            num = int(f[:-5].split('_')[-1])
            if nums:
                if num in nums:
                    snapfiles[num] = f
            else:
                snapfiles[num] = f
        return snapfiles

    def set_snapshots(self, *nums):
        self.snapfiles = self.find_snapshots(*nums)

    def load_snapshot(self, num, *load_keys,**kwargs):
        try:
            fname = self.snapfiles[num]
        except KeyError:
            try:
                fname = self.find_snapshots(num)[num]
            except KeyError:
                raise IOError('Sim ' + self.name + ' snapshot '
                              + str(num) + ' not found!')
        snap = snapshot.File(self, fname)

        if load_keys:
            snap.gas.load_data(*load_keys,**kwargs)
        return snap

    def multitask(self,plot_func,*data):
        maxprocs = mp.cpu_count()
        # Hack to take advantage of larger memory on r900 machines
        if 'r900' not in subprocess.check_output(['uname','-n']):
            maxprocs /= 2
        file_queue = Queue.Queue()
        data_queue = Queue.Queue(maxprocs/4)
        loader = snapshot.Loader(self.load_snapshot, file_queue, data_queue)
        loader.start()
        for snap in self.snapfiles.keys():
            args = (snap,)+data
            file_queue.put(args)
        file_queue.put(None)

        done = False
        while not done:
            jobs = []
            for i in xrange(maxprocs):
                snap = data_queue.get()
                if snap is None:
                    done = True
                    break
                else:
                    snap.close()
                    wp = self.plotpath + self.name
                    p = mp.Process(target=plot_func, args=(snap, wp))
                    jobs.append(p)
                    p.start()
            print '\nClearing Queue!\n'
            for process in jobs:
                process.join()
            print '\nQueue Cleared!\n'
        for process in jobs:
            process.join()

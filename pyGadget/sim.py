# sim.py
# Jacob Hummel
import os
import glob
import numpy
import Queue
import subprocess
import multiprocessing as mp

import sink
import units
import snapshot
import plotting

class Simulation(object):
    """
    Class for simulations.  Primarily for gathering metadata such as file
    paths and unit coversions.
    """
    def __init__(self, sim_name, **simargs):
        super(Simulation,self).__init__()
        self.name = sim_name
        self.datapath = simargs.pop('datapath',os.getenv('HOME')+'/sim/')
        self.plotpath = simargs.pop('plotpath',(os.getenv('HOME')
                                               +'/data/simplots/'))
        self.sinkpath = simargs.pop('sinkpath',os.getenv('HOME')+'/data/sinks/')
        self.hires = simargs.pop('refine', True)
        self.sink_tracking = simargs.pop('track_sinks', False)
        self.coordinates = simargs.pop('coordinates', 'physical')
        self.batch_viewscale = None

        self.units = units.Units(**simargs)
        self.units.set_coordinate_system(self.coordinates)

        self.snapfiles = self.find_snapshots()
        try:
            self.sink0 = sink.SingleSink(self.sinkpath + self.name)
            print "Found sinkfiles.  Loading sinkdata."
            self.tsink = self.sink0.time[0]
        except IOError:
            self.tsink = None

    def track_sinks(self, boolean=True):
        self.sink_tracking = boolean

    def refine_by_mass(self, boolean=True):
        self.hires = boolean

    def set_coordinate_system(self,coordinates):
        self.units.set_coordinate_system(coordinates)
        self.coordinates = self.units.coordinate_system

    def set_batch_viewscale(self, scale):
        unit = "".join(ch if not ch.isdigit() else "" for ch in scale)
        if unit in self.units.lengths.keys():
            self.batch_viewscale = scale
        else:
            raise KeyError

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
        if ((kwargs.pop('track_sinks',False)) or self.sink_tracking):
            kwargs['track_sinks'] = True
        if not self.hires:
            kwargs['refine_gas'] = False
        try:
            fname = self.snapfiles[num]
        except KeyError:
            try:
                fname = self.find_snapshots(num)[num]
            except KeyError:
                raise IOError('Sim ' + self.name + ' snapshot '
                              + str(num) + ' not found!')
        snap = snapshot.File(self, fname, **kwargs)

        if kwargs.get('track_sinks', False):
            if snap.gas._sink_indices.size > 0:
                if snap.sim.tsink is None:
                    snap.sim.tsink = snap.header.Time * units.Time_yr

        if load_keys:
            snap.gas.load_data(*load_keys,**kwargs)
        else:
            snap.gas.cleanup()
        return snap

    def multitask(self,plot_func,*data,**kwargs):
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
        if kwargs.pop('parallel', True):
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
        else:
            while not done:
                snap = data_queue.get()
                if snap is None:
                    done = True
                    break
                else:
                    snap.close()
                    wp = self.plotpath + self.name
                    plot_func(snap, wp)

#!/usr/bin/env python
# xdmf.py
# Jacob Hummel
"""
Script to generate xdmf wrappers for Gadget HDF5 files.
"""
import os
import glob
import subprocess

def hdfwrap(hdf5path=None, namestyle="snapshot_???.hdf5", 
            template='gadget.lxmf'):
    if not os.path.exists(template):
        raise IOError("Specified XDMF template file, '"+template+"' not found.")
    if not hdf5path:
        files = glob.glob(namestyle)
    else:
        if(hdf5path[-1] is not '/'): 
            if(namestyle[0] is not '/'):
                hdf5path = hdf5path + '/'
        files = glob.glob(hdf5path+namestyle)

    for hdf5file in files:
        xdmfile = open(hdf5file.replace('hdf5', 'xmf'), 'w')
        xdmfgen = subprocess.Popen(['XdmfGenerate', template, 
                                    hdf5file, '0'], stdout=xdmfile)
        xdmfile.close()

if __name__ == '__main__':
    hdfwrap()

#!/usr/bin/env python
# xdmf.py
# Jacob Hummel
"""
Script to generate xdmf wrappers for Gadget HDF5 files.
"""
import os
import glob
import subprocess
import pyGadget

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
    files.sort()

    outfiles = []
    for hdf5file in files:
        snap = pyGadget.snapshot.Load(hdf5file)
        proc = subprocess.Popen(['XdmfGenerate', template, 
                                 hdf5file, '0'], stdout=subprocess.PIPE)
        text = proc.communicate()[0]
        index = text.find(r'<Grid Name="PartType0" GridType="Uniform">')
        text = text[:index-2]+(r'<Time Type="Single" Value="'
                             +str(snap.header.Time)+'"/>\n')+text[index-6:]
        snap.close()

        xdmfile = open(hdf5file.replace('hdf5', 'xmf'), 'w')
        xdmfile.write(text)
        outfiles.append(xdmfile.name)
        xdmfile.close()
    return outfiles

def wrap_all(files):
    timexmf = open('series.xmf','w')
    timexmf.write('<?xml version="1.0" ?>\n')
    timexmf.write('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.1">\n')
    timexmf.write('  <Domain>\n')
    timexmf.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')
    for xdmfile in files:
        timexmf.write('    <xi:include href="'+xdmfile
                      +'" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n')
    timexmf.write('    </Grid>\n')
    timexmf.write('  </Domain>\n')
    timexmf.write('</Xdmf>')
    timexmf.close()
if __name__ == '__main__':
    files = hdfwrap()
    wrap_all(files)

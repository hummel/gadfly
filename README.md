GADFLY
========

Tools for reading and analyzing simulation snapshot files in the HDF5 format produced by the SPH simulation codes GADGET and GIZMO.

*A complete description of gadfly's capabilities can be found in this [paper](http://arxiv.org/abs/1603.05190), submitted to Publications of the Astronomical Society of the Pacific.*

Installation
------------
This code only runs on python 2, and is only tested on python 2.7. It requires the following packages:

  + numpy v1.7 or greater
  + pandas
  + h5py
  + **numba**
  + **scipy.weave**

The bolded entries are optional, needed only to take advantage of `gadfly`'s volume rendering capabilities, and only one of the two packages is needed. After cloning or otherwise downloading the repository, type 

`python setup.py install`

to install `gadfly`. 

Documentation
=============
Documentation is available on [Read the Docs](http://gadfly.readthedocs.org/en/latest/#), but is very much a work in progress at the moment.

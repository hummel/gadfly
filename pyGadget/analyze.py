# analyze.py
# Jacob Hummel
import numpy
import statistics
from . import units
from . import constants

#===============================================================================
def find_center(x, y, z, dens, **kwargs):
	centering = kwargs.pop('centering','avg')
	dens_limit = kwargs.pop('dens_limit', 1e11)
	nparticles = kwargs.pop('centering_npart', 100)
	verbose = kwargs.pop('centering_verbose', True)
	retcen = kwargs.pop('retcen',False)
	if centering == 'avg':
		hidens = numpy.where(dens >= dens_limit)[0]
		while hidens.size < nparticles:
		    dens_limit /= 2
		    hidens = numpy.where(dens >= dens_limit)[0]
		if verbose:
		    print ('Center averaged over %d particles' %nparticles)
		    print ('Center averaged over all particles with density '\
			       'greater than %.2e particles/cc' %dens_limit)
		#Center on highest density clump, rejecting outliers:
		cx = numpy.average(statistics.reject_outliers(x[hidens]))
		cy = numpy.average(statistics.reject_outliers(y[hidens]))
		cz = numpy.average(statistics.reject_outliers(z[hidens]))
		print 'Density averaged box center: %.3e %.3e %.3e' %(cx,cy,cz)
	elif centering == 'max':
		center = dens.argmax()
		cx,cy,cz = x[center], y[center], z[center]
		print 'Density maximum box center: %.3e %.3e %.3e' %(cx,cy,cz)
	else:
		cx = (x.max() + x.min())/2
		cy = (y.max() + y.min())/2
		cz = (z.max() + z.min())/2
		print 'Simple box center: %.3e %.3e %.3e' %(cx,cy,cz)
	x -= cx
	y -= cy
	z -= cz
        return x,y,z


# nbody.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 nbody (dark matter)
particle data.
"""
import hdf5

class PartTypeDM(hdf5.PartTypeX):
    """
    Class for Dark Matter particles.
    Available to extend class PartTypeX for Dark Matter specific applications.
    """
    def __init__(self, file_id, units, **kwargs):
        super(PartTypeDM,self).__init__(file_id,1, units)

        self._refined = None
        self.refine = kwargs.pop('refine_DM', False)
        if self.refine:
            print 'Turning on dark matter particle refinement.'
            self.locate_refined_particles()

    def locate_refined_particles(self):
        mass = self.get_masses()
        minimum = numpy.amin(mass)
        self._refined = numpy.where(mass <= minimum)[0]
        self.cleanup()

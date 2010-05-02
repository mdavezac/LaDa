""" Interface to c/fortran/c++ minimizers

    This is mostly for VFF. Should use numpy mininmizers where possible.
"""
from minimizer import *

def _print_minimizer(self):
  return "# Minimizer definition.\n"\
         "minimizer = Minimizer()\n"\
         "minimizer.type = \"%s\"\n"\
         "minimizer.tolerance = %e\n"\
         "minimizer.itermax = %i\n"\
         "minimizer.linetolerance = %e\n"\
         "minimizer.linestep = %e\n"\
         "minimizer.strategy = \"%s\"\n"\
         "minimizer.verbose = %s\n"\
         "minimizer.uncertainties = %e\n"\
         "minimizer.up = %i\n"\
         "minimizer.use_gradient = %s\n"\
         % ( self.type, self.tolerance, self.itermax, self.linetolerance, \
             self.linestep, self.strategy, "True" if self.verbose else "False", \
             self.uncertainties, self.up, "True" if self.use_gradient else "False" )
Minimizer.__str__ = _print_minimizer


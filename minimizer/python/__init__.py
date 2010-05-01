""" Interface to c/fortran/c++ minimizers

    This is mostly for VFF. Should use numpy mininmizers where possible.
"""
from minimizer import *

def _print_minimizer(self):
  return "# Minimizer definition.\n"\
         "minimizer = Minimizer()\n"\
         "minimizer.type = \"%s\"\n"\
         "minimizer.convergence = %e\n"\
         "minimizer.itermax = %i\n"\
         "minimizer.linetolerance = %e\n"\
         "minimizer.linestep = %e\n"\
         "minimizer.strategy = \"%s\"\n"\
         "minimizer.verbose = %s\n"\
         "minimizer.uncertainties = %e\n"\
         "minimizer.up = %i\n"\
         "minimizer.gradient = %s\n"\
         % ( minimizer.convergence, minimizer.itermax, minimizer.linetolerance, \
             minimizer.linestep, minimizer.strategy, "True" if minimizer.verbose else "False", \
             minimizer.uncertainties, minimizer.up, "True" if minimizer.gradient else "False" )
Minimizer.__str__ = _print_minimizer

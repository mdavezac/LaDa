""" Interface to c/fortran/c++ minimizers

    This is mostly for VFF. Should use numpy mininmizers where possible.
"""
__docformat__ = "restructuredtext en"
__all__ = ["Minimizer"]
      
class Minimizer(object):
  """ Interface to different C++ minimizers. """
  def __init__(self):
     """ Initializes minimizer wrapper. """
     self.type           = "gsl_bfgs2"
     self.tolerance      = 1e-6
     self.itermax        = 50
     self.linetolerance  = 1e-2
     self.linestep       = 1e-3
     self.strategy       = "fast"
     self.verbose        = False
     self.uncertainties  = 0.1
     self.up             = 1
     self.use_gradient   = True
     self.zeps           = 0.1

  @property 
  def type(self):
    """ The type of minimizer variant to use. """
    return self._type
  @type.setter
  def type(self, value):
    from _minimizer import all_types
    assert value.lower() in all_types, \
           RuntimeError("{0} is not a known/compiled minimizer. Use {1}.".format(name, all_types))
    self._type = value.lower()
  @property 
  def strategy(self):
    """ Strategy of MINUIT2 optimizers. """
    return self._strategy
  @strategy.setter
  def strategy(self, value):
    assert value.lower() in ["fast", "slow", "slowest"], \
           RuntimeError("strategy must be one of \"fast\", "
                        "\"slow\", \"slowest\" (vs {0}).".format(value))
    self._strategy = value.lower()

  @property
  def _minimizer(self):
    """ Creates a minimizer instance. """
    from _minimizer import Minimizer as ImplMinimizer
    result = ImplMinimizer()
    result.set( type = self.type,
                tolerance = self.tolerance,
                itermax = self.itermax,
                linetolerance = self.linetolerance,
                linestep = self.linestep,
                strategy = self.strategy,
                verbose = self.verbose, 
                uncertainties = self.zeps if self.type == "frprmn" else self.uncertainties,
                use_gradient = self.use_gradient, 
                up = self.up )
    return result

  def __call__(self, function, arg, **kwargs):
    """ Performs minimization.

        :Parameters:
          function 
            Callable bject (with `__call__` attribute) over which to perform minimization.
          arg 
            Starting points from where to minimize on input. Needs to be
            deepcopiable. It should not be an intrinsic type.
          kwargs 
            Attributes of the minimizer to modify here. 

        :return: A 2-tuple consisting of the minimized function value and the
                 locus of the minimum.
    """ 
    from copy import deepcopy
    try: result = deepcopy(arg)
    except: raise ValueError("arg should be deep-copiable")
    this = self.__class__()
    for key, value in kwargs:
      assert key[0] != '_', ValueError("Cannot modify private attribute " + key + ".")
      assert key in  this.__dict__, ValueError("Unknown attribute " + key + ".")
      setattribute(this, key, value)

    return this._minimizer(function, result), result

  def __repr__(self):
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
           "minimizer.zeps = %s\n"\
           "# End of minimizer definition.\n"\
           % ( self.type, self.tolerance, self.itermax, self.linetolerance, \
               self.linestep, self.strategy, "True" if self.verbose else "False", \
               self.uncertainties, self.up, "True" if self.use_gradient else "False",\
               self.zeps )
  def _copy_to_cpp(self, other):
    """ Copies minimizer to cpp type.

        The C++ instance is actually a boost::variant object. This means we
        have to jump through hoops to interface with it. 
    """ 
    self._minimizer._copy_to_cpp(other)


""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
from .extract import Extract as GAExtract

__all__ = ['Darwin', 'maximize', 'minimize']
def maximize(indivA, indivB, tolerance=1e-12):
  """ Compares two individuals. 

      This method can be used to *maximize* the objective function during the GA
      evolution.

      :Parameters:
        indivA 
          Individual with a scalar fitness attribute.
        indivB 
          Individual with a scalar fitness attribute.

      :return:
        - 1 if ``indivA.fitness < indivB.fitness``
        - 0 if ``indivA.fitness == indivB.fitness``
        - -1 if ``indivA.fitness > indivB.fitness``
  """
  from math import fabs
  a, b = indivA.fitness, indivB.fitness
  if hasattr(b, 'rescale') and hasattr(a, 'magnitude'):
    a, b = a.magnitude, b.rescale(a.units).magnitude
  if fabs(a-b) <  tolerance: return 0
  return 1 if float(a-b) < 0 else -1

def minimize(indivA, indivB, tolerance=1e-12):
  """ Compares two individuals. 

      This method can be used to *minimize* the objective function during the GA
      evolution.

      :Parameters:
        indivA 
          Individual with a scalar fitness attribute.
        indivB 
          Individual with a scalar fitness attribute.

      :return:
        - 1 if ``indivA.fitness > indivB.fitness``
        - 0 if ``indivA.fitness == indivB.fitness``
        - -1 if ``indivA.fitness < indivB.fitness``
  """
  from math import fabs
  a, b = indivA.fitness, indivB.fitness
  if hasattr(b, 'rescale') and hasattr(a, 'magnitude'):
    a, b = a.magnitude, b.rescale(a.units).magnitude
  if fabs(a-b) <  tolerance: return 0
  return -1 if float(a-b) < 0 else 1

class Darwin(object): 
  """ GA functional for optimizations of epitaxial structures. """
  Extract = GAExtract
  """ Holds all GA parameters """
  ordinals = GAExtract.ordinals
  """ Names of each historical age in the GA. """
  OUTCAR = GAExtract.OUTCAR
  """ Output filename. """
  ERRCAR = GAExtract.ERRCAR
  """ Error filename. """
  FUNCCAR = GAExtract.FUNCCAR
  """ Functional filename. """

  def __init__(self, evaluator, **kwargs): 
    """ Initializes a GA functional. 
         
        :Parameters:
          evaluator : `lada.ga.escan.elemental.Bandgap`
            Functional which uses vff and escan for evaluating physical properties.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
        :Kwarg pools:
          Number of pools over which to perform calculations. 
        :Kwarg rootworkdir:
          Root of the working directory where actual calculations are carried
          out. Defaults to $SCRACTH
        :Kwarg comparison:
          Function to compare the fitness of two individual. 
        :Kwarg history:
          If true, will use history in taboo.
        :KWarg histlimit:
          Number of invididuals kept in history. Limits the size of history to
          allow for faster processing.
    """
    from .history import History
    super(Darwin, self).__init__()

    self.evaluator = evaluator
    """ Evaluator object taking an individual and computing its fitness. """

    # other parameters
    self.popsize = kwargs.pop("popsize", 100)
    self.rate    = kwargs.pop("rate", max(1e0/(float(self.popsize)-0.5), 0.1))
    self.max_gen = kwargs.pop("max_generations", kwargs.pop("max_gen", 0))

    self.current_gen = kwargs.pop("current_gen", 0)
    """ Current generation. """

    self.population = []
    """ Current population """
    self.offspring = []
    """ Current offspring """

    self.age = None
    """ Current historical age (e.g. number of restarts). """

    self.pools = kwargs.pop("pools", 1)
    """ Number of pools of processors over which to perform calculations. """
    self.rootworkdir = kwargs.pop("rootworkdir", "$SCRACTH")
    """ Root directory where calculations are performed. """
    self.comparison = kwargs.pop("comparison", maximize)
    """ Comparison operator between individuals. """
    self.history = History() if kwargs.pop("history", False) else None 
    """ Holds list of visited candidates, if any. """
    if "histlimit" in kwargs: self.history.limit = kwargs.pop("histlimit")

  @property
  def rootworkdir(self):
    """ Root of the working directory where actual calculations are carried out. """
    return self.evaluator._outdir.envvar

  @rootworkdir.setter
  def rootworkdir(self, value):
    self.evaluator._outdir.envvar = value
    
  @property
  def checkpoints(self):
    """ Checkpoints functions.
    
        Defined as a property to avoid hassle with __getstate__
    """
    from .standard import print_offspring, average_fitness, best, print_population
    return [ print_offspring, 
             print_population,
             average_fitness,
             best,
             self.__class__.print_nb_evals,
             self.__class__.save,
             self.__class__.timing,
             self.__class__.end_of_time,
             self.__class__.flush_out ]

  def end_of_time(self): 
    """ Returns False if the maximum number of generations has been reached. """
    return self.current_gen < self.max_gen if self.max_gen >= 0 else True

  def flush_out(self):
    """ Tries to flush current output. """
    from sys import stdout
    from os import fsync
    stdout.flush()
    try: fsync(stdout)
    except: pass
    return True

  def evaluation(self):
    """ Evaluates population. """
    from .standard import population_evaluation
    population_evaluation(self, self.evaluator, self.pools, self.comm)

  def taboo(self, indiv):
    """ No two individuals in the population and the offspring are the same. """
    from itertools import chain
    if any(indiv == a for a in chain(self.population, self.offspring)): return True
    if getattr(self, "history", None) is not None: # check whether a history object exists.
      if self.history(indiv): return True
    return False

  def restart(self, outdir, comm = None):
    """ Saves current status. """
    from os.path import exists, join

    # sanity check.
    if self.age == self.ordinals[0] or self.age is None: return

    if comm is None and hasattr(self, "comm"): comm = self.comm
    
    # creates extraction object
    extract = self.Extract(outdir, comm)

    # checks for restart file.
    path = join(join(outdir, extract.current_age), self.FUNCCAR)
    if not exists(path): 
      error = "Could not find restart file {0}.\n"\
              "Yet directory {1} exits.\nAborting\n"\
              .format(path, join(outdir, name))
      if self.do_print: print error
      raise RuntimeError(error)

    # copies populations and friends.
    self.population  = extract.functional.population
    self.offspring   = extract.functional.offspring
    self.current_gen = extract.functional.current_gen

    if self.do_print:
      print "Restarting from file {0}.\n"\
            "Restarting at generation {1}.\n"\
            .format(self.FUNCCAR, self.current_gen)
      if len(self.population) == 0: print "Restarting with new population.\n"
      else:
        print "Restarting with population:"
        for indiv in self.population: print indiv, indiv.fitness

      if len(self.offspring) == 0: print "Restarting with empty offspring.\n"
      else:
        print "Restarting with offspring:"
        for indiv in self.offspring: print indiv
      print "\nStarting {0} GA age.\n".format(self.age)

  def save(self):
    """ Saves current status. """
    from .standard import append_population
    from pickle import dump

    # only one proc should print.
    is_root = self.comm.rank == 0 if hasattr(self, "comm") else True
    if is_root:
      with open(self.FUNCCAR, "wb") as file: dump(self, file)
    append_population(self, self.population, self.Extract.OFFCAR)
    return True

  @property
  def do_print(self):
    """ Wether this process prints. """
    if not hasattr(self, "comm"): return True
    if self.comm is None: return True
    return self.comm.rank == 0

    
  @property
  def color(self):
    """ Returns color of this process and None if not pooled or MPI. """
    if self.pools < 2:            return None
    if not hasattr(self, "comm"): return None
    if self.comm is None:         return None
    pools = self.pools if self.comm.size >= self.pools else self.comm.size
    return self.comm.rank % pools


  def print_nb_evals(self):
    """ Prints current number of evaluations. """
    if self.color is not None:
      from ..mpi import all_reduce
      local_comm = self.comm.split(self.color)
      heads_comm = self.comm.split(1 if local_comm.rank == 0 else 2)
      nbcalc = all_reduce(heads_comm, getattr(self.evaluator, "nbcalc", 0), lambda x,y: x+y)
    else: nbcalc = getattr(self.evaluator, "nbcalc", 0)
    
    if self.do_print: print "  Number of functional evaluations: ", nbcalc

  def timing(self):
    """ Prints timing at each generation. """
    import time
    if not self.do_print: return True
    t = time.time() - self.start_time
    hour = int(float(t/3600e0))
    minute = int(float((t - hour*3600)/60e0))
    second = (t - hour*3600-minute*60)
    print "  Elapsed time: {0}:{1}:{2:.4f}.".format(hour, minute, second)
    return True


  def __call__(self, comm = None, outdir = None, inplace=False, **kwargs):
    """ Runs/Restart  GA 

        :Parameters:
          comm : mpi.Communicator
            MPI communicator when doing mpi. None otherwise.
          outdir : str
            Path to the output directory. Defaults to the current working
            directory.
          inplace : bool
            Ignored. Where calculations are performed depends on rootworkdir. 
    """
    import time
    from os import getcwd
    from os.path import join
    from copy import deepcopy
    from ..mpi import world
    from . import darwin as search
    from ..opt import redirect, RelativeDirectory, Changedir

    # make call stateless.
    this = deepcopy(self)

    local_time = time.localtime() 
    this.start_time = time.time() 
    if outdir is None: outdir = getcwd()
    outdir = RelativeDirectory(outdir, envvar=getcwd())
    # mpi stuff
    this.comm = comm if comm is not None else world
    this.comm.do_print = this.do_print

    # takes care of keyword arguments:
    kwargs.pop("external", None)
    if len(kwargs.keys()) > 0: 
      for key, value in kwargs.items():
        if hasattr(this, key): setattr(this, key, value)
        elif hasattr(this.evaluator, key): setattr(this.evaluator, key, value)
        else: assert hasattr(this, key), TypeError("Unknown argument {0}.".format(key))

    # gets current age
    this.age = this.Extract(outdir.path, comm).next_age
    # sets directory for calculations according to newly read age.
    if hasattr(this.evaluator, "outdir"): 
      this.evaluator._outdir.envvar = self.rootworkdir
      this.evaluator._outdir.relative = outdir.relative
      this.evaluator.outdir = join(this.evaluator.outdir, this.age)
    # sets directory for history file. This should be conserved across runs,
    # hence it is not on the scratch.
    if this.history is not None:
      this.history.directory = outdir.path 
      this.history.remove_stale(comm)


    # now goes to work
    with Changedir(join(outdir.path, this.age), comm=comm) as cwd:
      pyout = this.OUTCAR if this.do_print else '/dev/null' # this.OUTCAR + str(comm.rank) # 
      pyerr = this.ERRCAR if this.do_print else '/dev/null' # this.OUTCAR + str(comm.rank) # 
      with redirect(pyout=pyout, pyerr=pyerr) as streams:
        if this.do_print:
          print "# GA calculation on ", time.strftime("%m/%d/%y", local_time),\
                " at ", time.strftime("%I:%M:%S %p", local_time)
        # reloads if necessary
        if this.age != this.ordinals[0]: this.restart(outdir.path, comm=comm)
        # call back for further setting up in derived class.
        this._further_setup()
        # runs.
        search.run(this)
        this.timing()
        
  def _further_setup(self): 
    """ Derived class can override this function for last minute setup. """
    pass

  def __str__(self):
    return "Replacement rate: %f\nPopulation size: %i\nMaximum number of generations: %i\n" \
           % (self.rate, self.popsize, self.max_gen)

  def __getstate__(self):
    """ Returns current state. """
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d

  def __setstate__(self, state):
    """ Resets current state. """
    self.__dict__.update(state)
    # takes care of older pickles.
    if not hasattr(self, "comparison"): self.comparison = maximize
    if not hasattr(self, "history"):    self.history    = None

  def __repr__(self):
    """ Returns representation of this instance. """
    return "from {0} import {1}\n{2}\n"\
           "from {3.comparison.__module__} import {3.comparison.__name__}\n"\
           "from {3.history.__module__} import {3.history.__name__}\n"\
           "functional = {1}(evaluator)\n"\
           "functional.pools       = {3.pools}\n"\
           "functional.cm_rate     = {3.cm_rate}\n"\
           "functional.popsize     = {3.popsize}\n"\
           "functional.max_gen     = {3.max_gen}\n"\
           "functional.age         = {4}\n"\
           "functional.current_gen = {3.current_gen}\n"\
           "functional.comparison  = {3.comparison.__name__}\n"\
           "functional.history     = {5}\n"\
           "# functional.rootworkdir, please set evaluator.outdir instead\n"\
           .format( self.__class__.__module__, 
                    self.__class__.__name__,
                    repr(self.evaluator),
                    self, repr(self.age), repr(self.history) )




        




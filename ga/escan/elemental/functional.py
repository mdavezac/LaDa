""" GA as a functional, for use with ladajobs. """
__docformat__ = "restructuredtext en"
from .extract import Extract as GAExtract

__all__ = ['Darwin']

class Darwin: 
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
        :Kwarg crossover:
          Crossover object. Defaults to `lada.ga.escan.elemental.Crossover` ``(rate=0.5)``.
        :Kwarg mutation:
          Mutation object. `lada.ga.elemental.Mutation` ``(rate=2e0/float(len(.evaluator))``
        :Kwarg cm_rate:
          Crossover vs Mutation rate. Defaults to 0.1.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg mean_conc:
          Mean concentration of individuals on random initialization. 
        :Kwarg stddev_conc:
          Standard deviation of concentration of individuals on random initialization.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
    """
    from .. import CompareSuperCells
    from . import Crossover, Mutation
    from ....opt import RelativeDirectory


    self.evaluator = evaluator
    """ Evaluator object taking an individual and computing its fitness. """
    lattice = self.evaluator.escan.vff.lattice
    unit_cell = self.evaluator.converter.structure.cell
    self.compare = CompareSuperCells(lattice, unit_cell, self.evaluator.converter)
    """ Compares two different individuals. """


    # mating operations
    self.cm_rate = kwargs.pop("cm_rate", 0.1)
    """ Crossover versus Mutation rate. """
    self.crossover = kwargs.pop("crossover", Crossover(rate=0.5))
    """ Crossover operator. """
    self.mutation = kwargs.pop("mutation", Mutation(rate=2e0/float(len(self.evaluator))))
    """ Mutation operator. """

    # other parameters
    self.popsize = kwargs.pop("popsize", 100)
    self.rate    = kwargs.pop("rate", max(1e0/(float(self.popsize)-0.5), 0.1))
    self.max_gen = kwargs.pop("max_generations", kwargs.pop("max_gen", 0))

    # Parameters for starting population.
    self.mean_conc = kwargs.pop("mean_conc", 0.5)
    """ Mean concentration of individuals on random initialization. """
    self.stddev_conc = kwargs.pop("stddev_conc", 0.5)
    """ Standard deviation of concentration of individuals on random initialization. """

    self.current_gen = kwargs.pop("current_gen", 0)
    """ Current generation. """

    self.population = []
    """ Current population """
    self.offspring = []
    """ Current offspring """

    self.age = None
    """ Current historical age (e.g. number of restarts). """
    
  @property
  def checkpoints(self):
    """ Checkpoints functions.
    
        Defined as a property to avoid hassle with __getstate__
    """
    from ...standard import print_offspring, average_fitness, best
    return [ print_offspring, 
             average_fitness,
             best,
             self.__class__.print_nb_evals,
             self.__class__.save,
             self.__class__.timing,
             self.__class__.check_generations,
             self.__class__.flush_out ]

  def check_generations(self): 
    """ Returns False if the maximum number of generations has been reached. """
    return self.current_gen < self.max_gen if self.max_gen >= 0 else True

  def flush_out(self):
    """ Tries to flush current output. """
    from sys import stdout
    from os import fsync
    stdout.flush()
    try: fsync(stdout)
    except: pass

  def evaluation(self):
    """ Evaluates population. """
    from ...standard import population_evaluation
    population_evaluation(self, self.evaluator, self.comm, self.pools)

  def mating(self, *args, **kwargs):
    """ Calls mating operations. Removes hash code. """
    from .. import CompareSuperCells
    from ...standard import Mating
     
    # construct sequential mating operator.
    mating = Mating(sequential=False)
    mating.add(self.crossover, rate=self.cm_rate)
    mating.add(self.mutation,  rate=1-self.cm_rate) 
    # calls mating operator.
    result = mating(self, *args, **kwargs)
    # removes any has-code.
    CompareSuperCells.remove_hashcode(result)
    return result


  def taboo(self, indiv):
    """ No two individuals in the population and the offspring are the same. """
    from itertools import chain
    for a in chain(self.population, self.offspring):
      if self.compare(a, indiv): return True
    return False

  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    from numpy.random import normal, shuffle
    result = Individual(size=len(self.evaluator))
    if self.mean_conc == None or self.stddev_conc == None: return result
    N = 0
    while N == result.genes.size or N == 0:
      x = normal(self.mean_conc, self.stddev_conc)
      if x >= 1 or x <= 0e0: continue
      N = int(x * float(result.genes.size))
    result.genes[:N] = 1 # second type of appropriate site in lattice.
    result.genes[N:] = 0 # first type of appropriate site in lattice.
    shuffle(result.genes)
    return result

  @staticmethod
  def cmp_indiv(a, b, tolerance = 1e-12 ):
    """ Compares two individuals. 

        Maximizes fitness by default.
    """
    from math import fabs
    if fabs(a.fitness - b.fitness) <  tolerance: return 0
    return 1 if a.fitness < b.fitness else -1

  def restart(self, outdir, comm = None):
    """ Saves current status. """
    import cPickle
    from os.path import exists, join

    # sanity check.
    if self.age == self.ordinals[0] or self.age == None: return

    if comm == None and hasattr(self, "comm"): comm = self.comm
    
    # creates extraction object
    extract = self.Extract(outdir, comm)

    # checks for restart file.
    path = join(join(outdir, extract.current_age), self.FUNCCAR)
    if not exists(path): 
      error = "Could not find restart file {0}.\n"\
              "Yet directory {1} exits.\nAborting\n"\
              .format(path, join(outdir, name))
      if self.do_print: print error
      raise RuntimeErorr(error)

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
    from pickle import dump

    # only one proc should print.
    is_root = self.comm.rank == 0 if hasattr(self, "comm") else True
    if is_root:
      with open(self.FUNCCAR, "wb") as file: dump(self, file)
    return True

  @property
  def do_print(self):
    """ Wether this process prints. """
    if not hasattr(self, "comm"): return True
    if self.comm == None: return True
    return self.comm.rank == 0

    
  @property
  def color(self):
    """ Returns color of this process and None if not pooled or MPI. """
    if self.pools < 2:            return None
    if not hasattr(self, "comm"): return None
    if self.comm == None:         return None
    pools = self.pools if self.comm.size >= self.pools else self.comm.size
    return self.comm.rank % pools


  def print_nb_evals(self):
    """ Prints current number of evaluations. """
    is_mpi = self.comm != None
    do_print = 0 if not is_mpi else self.comm.do_print 
    if self.color != None:
      from boost.mpi import all_reduce
      local_comm = self.comm.split(self.color)
      heads_comm = self.comm.split(1 if local_comm.rank == 0 else 2)
      nbcalc = all_reduce(heads_comm, self.evaluator.nbcalc, lambda x,y: x+y)
    else: nbcalc = self.evaluator.nbcalc
    
    if self.do_print:
      print "  Number of functional evaluations: ", nbcalc

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


  def __call__(self, comm = None, outdir = None, inplace=None, **kwargs):
    """ Runs/Restart  GA """
    import time
    from os import getcwd
    from os.path import join
    from copy import deepcopy
    from boost.mpi import world
    from ... import darwin as search
    from ....opt import redirect, RelativeDirectory, Changedir

    local_time = time.localtime() 
    self.start_time = time.time() 
    outdir = RelativeDirectory(outdir)
    if outdir == None: outdir = getcwd()
    # mpi stuff
    self.comm = comm if comm != None else world
    self.comm.do_print = self.do_print

    # takes care of keyword arguments:
    if len(kwargs.keys()) > 0: 
      this = deepcopy(self)
      for key, value in kwargs.items():
        assert hasattr(self, key), TypeError("Unknown argument {0}.".format(key))
        setattr(self, key, value)
      this(comm=comm)
      return

    # gets current age
    self.age = self.Extract(outdir.path, comm).next_age
    self.evaluator._outdir.relative = outdir.relative
    self.evaluator.outdir = join(self.evaluator.outdir, self.age)
#   if self.color != None: 
#     self.evaluator.outdir = join(self.evaluator.outdir, "pool_{0}".format(self.color))

    # now goes to work
    with Changedir(join(outdir.path, self.age), comm=comm) as cwd:
      pyout = self.OUTCAR if self.do_print else '/dev/null'
      pyerr = self.ERRCAR if self.do_print else '/dev/null'
      with redirect(pyout=pyout, pyerr=pyerr) as streams:
        if self.do_print:
          print "# GA calculation on ", time.strftime("%m/%d/%y", local_time),\
                " at ", time.strftime("%I:%M:%S %p", local_time)
        # reloads if necessary
        if self.age != self.ordinals[0]: self.restart(outdir.path, comm=comm)
        # runs.
        search.run(self)
        self.timing()
        
    del self.comm

  def __str__(self):
    return "Replacement rate: %f\nPopulation size: %i\nMaximum number of generations: %i\n" \
           "Size of each individuals: %i\n"\
           % (self.rate, self.popsize, self.max_gen, len(self.evaluator))

  def __getstate__(self):
    """ Returns current state. """
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d

  def __setstate__(self, state):
    """ Resets current state. """
    self.__dict__.update(state)

  def __repr__(self):
    """ Returns representation of this instance. """
    return "from {0} import {1}\n{2}\n{3}\n{4}\n\n"\
           "functional = {1}(evaluator, crossover=crossover, mutation=mutation)\n"\
           "functional.pools       = {5.pools}\n"\
           "functional.cm_rate     = {5.cm_rate}\n"\
           "functional.popsize     = {5.popsize}\n"\
           "functional.max_gen     = {5.max_gen}\n"\
           "functional.age         = {5.age}\n"\
           "functional.current_gen = {5.current_gen}\n"\
           "functional.mean_conc   = {5.mean_conc}\n"\
           "functional.stddev_conc = {5.stddev_conc}\n"\
           .format( self.__class__.__module__, 
                    self.__class__.__name__,
                    repr(self.evaluator),
                    repr(self.crossover).replace("ga_operator", "crossover"),
                    repr(self.mutation).replace("ga_operator", "mutation"),
                    self )




        




""" Holds standard genetic algorithm operations. 

    :group Checkpoints: best, print_population, print_offspring,
                        average_fitness, _check_generation
"""
__docformat__ = "restructuredtext en"
def bound_method(self, method):
  """ Returns a method bound to self. """
  from new import instancemethod
  return instancemethod(method, self, self.__class__) 

class Taboo(object):
  """ A container of taboo operators.

      By default, a diversity taboo operator is added. 
  """

  def __init__(self, diversity=True):
    """ Creates a Taboo container.

        :Param diversity: if True then a diversity constraint is added to the
          Taboo container. If a callable, uses should be a comparison operator
          to use with the diversity operator.  Otherwise the container is empty
          on initialization.
        :type diversity: Boolean, or callable. 
    """
    super(Taboo, self).__init__()
    self.taboos = []


    if hasattr(diversity, "__call__"):
      
      def diversity_taboo(self, indiv):
        """ taboo makes sure that no two individuals in the population and the
            offspring are the same. """
        from itertools import chain
        comparison_operator = diversity
        for a in chain(self.population, self.offspring):
          if comparison_operator(a, indiv): return True
        return False

      self.taboos.append(diversity_taboo)

    elif diversity:

      def diversity_taboo(self, indiv):
        """ taboo makes sure that no two individuals in the population and the
            offspring are the same. """
        return indiv in self.population or indiv in self.offspring

      self.taboos.append(diversity_taboo)


  def add(self, taboo):
    """ Adds a taboo operator to the list.

        A taboo operator takes the darwin class and the individual as arguments
        and returns True if the individual is taboo.
    """
    self.taboos.append(taboo)

  def __call__(self, darwin, indiv):
    """ Returns true if any one operator returns true. """
    import sys
    for taboo in self.taboos:
      if taboo(darwin, indiv): return True

    return False

def tournament( self, size = 2 ):
  """ deterministic tournament """
  import random
  list_ = range(len(self.population))
  random.shuffle(list_)
  list_ = list_[:size]
  result = 0
  for b in list_:
    if self.cmp_indiv(self.population[b], self.population[result]):
      result = b;
  return result

def cmp_indiv( a, b, tolerance = 1e-12 ):
  """ Compares two individuals. 

      Minimizes fitness by default.
  """
  from math import fabs
  if fabs(a.fitness - b.fitness) <  tolerance: return 0
  elif a.fitness < b.fitness: return -1
  return 1

def average_fitness(self):
  """ Prints out average fitness. """
  if not self.comm.do_print: return True;
  result = 0e0
  for indiv in self.population:
    result += indiv.fitness
  print "  Average Fitness: ", result / float(len(self.population))
  return True

def best(self):
  """ Checkpoint which prints out the best individual. """
  if not self.comm.do_print: return True
  best = None
  for indiv in self.population:
    if best == None or  self.cmp_indiv( best, indiv ) == 1: 
      best = indiv
  print "  Best Individual: ", best, best.fitness
  return True

def print_population(self):
  if not self.comm.do_print: return True
  print "  Population: "
  for indiv in self.population:
    print "    ", indiv, indiv.fitness
  return True

def print_offspring(self):
  if not self.comm.do_print: return True
  print "  Offspring: "
  for indiv in self.population:
    if indiv.birth == self.current_gen - 1: 
      print "    ", indiv, indiv.fitness
  return True

def _check_generation( self ):
  """ returns false if maximum number of generations was passed. 
      
      :attention: This checkpoint is always added by default. Users need not include it.
  """
  if self.max_gen < 0: return True
  return self.current_gen < self.max_gen
  
def serial_population_evaluation(self, evaluator, comm = None):
  """ Evaluates population and offspring serially.
  
      :Param evaluator: 
        Functional which performs actual calculations.

      Only individuals without a ``fitness`` attribute are evaluated. 
  """
  for indiv in self.population:
    if not hasattr(indiv, "fitness" ): 
      indiv.fitness = evaluator(indiv, comm=comm)
  for indiv in self.offspring:
    if not hasattr(indiv, "fitness" ): 
      indiv.fitness = evaluator(indiv, comm=comm)

def bleeder_evaluation(self, evaluator, pools, comm):
  """ MPI Evaluation using an on-disk dictionary.
  
      :Parameters: 
        self 
          A *darwin* type object with ``offspring`` and ``population``
          attributes. Only those on the root process matter. However, at the
          end of the call, all processes will share the same ``offspring`` and
          ``population``. 
        evaluator 
          An object capable of evaluating individuals.
        pools
          The number of pools of processes to create.
        comm  : None or boost.mpi.communicator
          Communication object. 
  """
  from operator import itemgetter
  from ..jobs import JobDict, Bleeder
  from ..mpi import Communicator

  comm = Communicator(comm)
  # creates job-dictionary.
  jobdict = JobDict()
    
  # goes throught individuals which need be evaluated
  # bleeder only looks at jobdict from root processs.
  if comm.is_root:
    for name, pop in [('off', self.offspring), ('pop', self.population)]:
      for index, indiv in enumerate(pop):
        job = jobdict / '{0}/{1}'.format(name, index)
        job.functional = evaluator
        job.jobparams['indiv'] = indiv
        job.jobparams['overwrite'] = True
        # lets not recompute already known individuals. Need to keep them as record though.
        if hasattr(indiv, "fitness"): job.tag() 

  # Creates a bleeding job-dict, eg a JobDictionary on file which can be
  # iterated over once, with each job being seen by only one pool of processes.
  # Modified jobs are saved by the bleeder!
  bleeder = Bleeder(jobdict, pools, comm)
  for job in bleeder:
    fitness = evaluator(job.indiv, outdir=job.name[1:], comm=bleeder.local_comm)
    job.indiv.fitness = bleeder.local_comm.broadcast(fitness)
  
  # recovers population and offspring.
  jobdict = bleeder.cleanup()
  if 'off' in jobdict: self.offspring  = [u.indiv for u in jobdict['off'].values()]
  if 'pop' in jobdict: self.population = [u.indiv for u in jobdict['pop'].values()]
    

def mpi_population_evaluation(self, evaluator, pools, comm = None):
  """ MPI Population and offspring evaluation.
  
      :Param evaluator: 
        Functional which performs actual calculations.
      :Param pools:
        Number of pools of processors across which to split evaluations.
      :Param comm: 
        group boost.mpi.communicator.

      Only individuals without a ``fitness`` attribute are evaluated. 
  """
  from boost.mpi import broadcast, scatter, all_gather, world
  from itertools import chain
  # split communicator along number of pools
  if pools == None: pools = comm.size
  if pools > comm.size: pools = comm.size
  color = comm.rank % pools
  local_comm = comm.split(color)
  heads_comm = comm.split(1 if local_comm.rank == 0 else 2)

  def check_pops(this, population):
    if not __debug__: return
    new_pop = broadcast(this.comm, population, 0)
    assert len(new_pop) == len(population),\
           RuntimeError("Populations across processes have different lengths.")
    for a, b in zip(new_pop, population):
      assert a == b, RuntimeError("Populations are not equivalent across processes.")
      ahas, bhas = hasattr(a, "fitness"), hasattr(b, "fitness")
      assert (ahas and bhas) or not (ahas or bhas),\
            RuntimeError("Populations do not have equivalently evaluated invidiuals")
      if ahas and bhas:
        assert a.fitness == b.fitness,\
        RuntimeError("Fitness are not equivalent across processes.")
    this.comm.barrier()


  gather_these = []
  # Now goes throught individuals which need be evaluated
  for index, indiv in enumerate(chain(self.population, self.offspring)):
    if hasattr(indiv, "fitness"): continue
    if index % pools == color: 
      fitness = evaluator(indiv, comm = local_comm)
      if local_comm.is_root: gather_these.append( (indiv, fitness) )

  # gathers all newly computed individuals. 
  if local_comm.is_root: gather_these = heads_comm.all_gather(gather_these)
  gather_these = local_comm.broadcast(gather_these)

  # now reinserts them into populations.
  for index, indiv in enumerate(chain(self.population, self.offspring)):
    if hasattr(indiv, "fitness"): continue
    assert len(gather_these) > index % pools, \
           RuntimeError("%s > %i %% %i" % (len(gather_these), index, pools))
    assert len(gather_these[index % pools]) != 0
    a, fitness = gather_these[index % pools].pop(0)
    indiv.__dict__.update(a.__dict__)
    indiv.fitness = fitness
  comm.barrier()
  for index, indiv in enumerate(chain(self.population, self.offspring)):
    if hasattr(indiv, "fitness"): continue
    assert False, "should not be here"
  comm.barrier()

  check_pops(self, self.population)
  check_pops(self, self.offspring)
  world.barrier()

def population_evaluation(self, evaluator, comm=None, pools=None):
  """ Chooses between MPI and serial evaluation. """
  from ..mpi import Communicator
  is_serial = pools == 1 or (comm.size == 1 if comm != None else True)
  comm = Communicator(comm)
  if is_serial: serial_population_evaluation(self, evaluator, comm = comm)
  else:         bleeder_evaluation(self, evaluator, min(pools, comm.size), comm)


class Mating(object):
  """ Aggregator of mating operators. 

      Mating operations can be added using the Mating.add bound method.
      Operations are called either sequentially[==and] or proportionally[==or] (either all
      or only one can be called). 

      Mating operations can be unary, binary, or ternary. The first argument passed
      to each operation is the new offspring, passed by reference. The other
      arguments are references to individual in the population. They should be
      treated as read-only.
  """
  def __init__(self, sequential=False): 
    """ Initializes the container of mating operator. """
    self.operators = []
    self.sequential = sequential
    """ If true, operators are applied sequentially. """

  def __repr__(self):
    string =  "from {0} import {1}".format(self.__class__.__module__, self.__class__.__name__)
    string +=  "mating = {0}({1})".format(self.__class__.__name__, self.sequential)
    for function, rate, n in self.operators:
      string += "mating.add({0}, {1})".format(repr(self.function), self.rate)
    return string


  def add(self, function, rate=1):
    """ Adds a mating operator, with a given rate, to the current list. """
    # checks for number of arguments to function.
    nb_args = -1
    if hasattr(function, "__class__") and issubclass(function.__class__, Mating): pass
    elif hasattr(function, "func_code"): nb_args = function.func_code.co_argcount
    else: nb_args = function.__call__.im_func.func_code.co_argcount - 1

    assert rate > 0e0, ValueError("rate argument cannot be negative (%s)." % (rate))
    self.operators.append( (function, rate, nb_args) )
  

  def __call__(self, darwin):
    """ Creates an offspring. """
    import random

    def call_function(function, n, indiv = None):
      """ Calls functions/functors mating operations. """
      from copy import deepcopy
      # calls other Mating instances.
      if n == -1: return function(darwin)

      individuals = []
      if indiv != None: individuals.append(indiv)
      else: individuals.append( deepcopy(darwin.population[darwin.selection(darwin)]) )

      # calls unaries
      if   n == 1: return function( individuals[0] )

      # calls binaries
      b = individuals[0]
      while( b not in individuals ): b = darwin.population[darwin.selection(darwin)]
      individuals.append(b)
      if n == 2: return function( individuals[0], individuals[1] )

      # calls ternaries
      b = individuals[0]
      while( b not in individuals ): b = darwin.population[darwin.selection(darwin)]
      individuals.append(b)
      if n == 3: return function( individuals[0], individuals[1], individuals[2] )

      raise "Mating operations has to be unary, binary, or ternary.\n"

    indiv = None
    if self.sequential: # choose any operator depending on rate.
      while indiv == None: # makes sure we don't bypass all mating operations
        for function, rate, n in self.operators:
          if random.random() < rate: indiv = call_function( function, n, indiv )
    else: # choose only one operator.
      max = 0e0
      for function, rate, n in self.operators: max += rate
      assert rate > 0e0;
      last = 0e0
      r = random.random() * max
      for function, rate, n in self.operators:
        if r <= last + rate:
          indiv = call_function( function, n )
          break
        last += rate

    assert indiv != None, "%s" % (self.sequential)
    return indiv


def add_checkpoint(self, _chk):
  """ Adds a checkpoint to self.checkpoints. """
  try: self.checkpoints.append( _chk ) 
  except AttributeError: self.checkpoints = [_chk]
  return self;

def fill_attributes(self):
  """ Checks self for correct attributes.

      Fills in where possible:

      - "taboo" defaults to standard.taboo
      - "selection" defaults to standard.selection
      - "population" defaults to empty list []
      - "popsize" defaults to len(self.population)
      - "offspring" defaults to empty list []
      - "rate" defaults to len(self.offspring)/self.popsize
      - standard._check_generation is ALWAYS added to the checkpoints
      - "current_gen" defaults to 0 
      - "max_gen" defaults to 100, or current_gen+100 if max_gen > current_gen.
      - "cmp_indiv" defaults to standard.cmp_indiv
  """
  import darwin 
  # must have an evaluation function.
  assert hasattr(self, "evaluation"), "No evaluation function!" 

  # Checks that self has an object Individual
  if not hasattr(self, "Individual"):
    from bitstring import Individual
    self.Individual = Individual

  # Checks whether self has a taboo object.
  if not hasattr(self, "taboo"): self.taboo = bound_method(self, Taboo())

  # Checks whether self has a selection object.
  if not hasattr(self, "selection"): self.selection = tournament

  # checks whether there is a population.
  if not hasattr(self, "population"): self.population = []
  
  # checks whether there is a popsize.
  if not hasattr(self, "popsize"): self.popsize = len(self.population)

  # makes sure we have something to do.
  assert self.popsize != 0, "population or popsize attributes required on input."

  # checks whether there is a population.
  if not hasattr(self, "offspring"): self.offspring = []
  
  # checks whether there is a popsize.
  if not hasattr(self, "rate"):
    self.rate = float(len(self.offspring)) / float(self.popsize)

  # makes sure we have something to do.
  assert self.rate > float(0), "offspring or rate attributes required on input."
 
  # checks whether there is a checkpoint.
  self = add_checkpoint(self, _check_generation)

  # checks current generation.
  if not hasattr(self, "current_gen"): self.current_gen = 0
  elif not hasattr(self, "max_gen"): self.max_gen = 100
  elif self.max_gen not in [0, -1] and self.max_gen < self.current_gen: self.max_gen += 100

  if not hasattr(self, "max_gen"): self.max_gen = self.current_gen + 100

  if not hasattr(self, "cmp_indiv"): self.cmp_indiv = cmp_indiv
  return self

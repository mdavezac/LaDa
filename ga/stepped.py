from abc import ABCMeta, abstractmethod

class Extract(object):
  """ Extraction class for StepGA. """
  def __init__(self, directory=None):
    """ Initializes an extraction object. 
  
        :param directory: 
          Directory where the SHELVECAR can be found.
          Defaults to current working directory.

        :raises RuntimeError: if no appropriate SHELVERCAR can be found.
    """
    from os.path import exists, isfile, isdir
    from os import makedirs
    from shelve import open as shelve_open
    from ..opt import RelativeDirectory
    super(Extract, self).__init__()

    self._directory = RelativeDirectory(directory)
    if not exists(self.directory):
      raise RuntimeError('{0} could not be found.'.format(self.directory))
    if not isdir(self.directory):
      raise RuntimeError('{0} is not a directory.'.format(self.directory))
    if not exists(self.shelvepath): 
      raise RuntimeError('{0} could not be found.'.format(self.shelvepath))
    if not isfile(self.shelvepath): 
      raise RuntimeError('{0} is not a file.'.format(self.shelvepath))

    try:
      if not exists(self.directory): makedirs(self.directory)
      shelve = shelve_open(self.shelvepath)
      if set(shelve.keys()) !=  set(['individuals', 'functionals', 'removed', 'added', 'new']):
        raise RuntimeError('{0} is not a StepGA SHELVECAR file.'.format(self.shelvepath))
    finally: shelve.close()

  @property
  def directory(self):
    """ Directory where output should be found. """
    return self._directory.path
  @property
  def shelvepath(self):
    """ Directory where output should be found. """
    from os.path import join
    return join(self._directory.path, StepGA.SHELVECAR)

  @property
  def functionals(self):
    """ Functionals at each generation. """
    from shelve import open as shelve_open
    try:
      shelve = shelve_open(self.shelvepath)
      return shelve['functionals']
    finally: shelve.close()

  @property 
  def populations(self):
    """ Populations at each generation. """
    from shelve import open as shelve_open
    try:
      shelve = shelve_open(self.shelvepath)
      return [func.population for func in self.functionals]
    finally: shelve.close()

  @property 
  def added(self):
    """ Individuals added at each generation. """
    from shelve import open as shelve_open
    try:
      shelve = shelve_open(self.shelvepath)
      result = []
      for add in shelve['added']:
        result.extend(shelve['individuals'][i] for i in add)
      return result
    finally: shelve.close()

  @property 
  def removed(self):
    """ Individuals removed at each generation. """
    from shelve import open as shelve_open
    try:
      shelve = shelve_open(self.shelvepath)
      result = []
      for r in shelve['removed']:
        result.extend(shelve['individuals'][i] for i in r)
      return result
    finally: shelve.close()

  def best(self, n=10): 
    """ Current n best individuals. 

        :param int n: 
           Number of individuals to print.
    """
    from operator import itemgetter
    from ..jobs.forwarding_dict import ForwardingDict

    individuals = self.individuals

    fitness = sorted( [(name, u.fitness) for name, u in individuals.iteritems()],\
                      key=itemgetter(1) )[:n]
    result = ForwardingDict(readonly=True, ordered=False)
    for key, value in fitness: result[key] = individuals[key]
    return result
    
  @property
  def individuals(self):
    """ Dictionary containing all individuals. """
    from shelve import open as shelve_open
    try: shelve = shelve_open(self.shelvepath)
    except: raise
    else: return shelve['individuals']
    finally: shelve.close()



class StepGA(object):
  """ Genetic Algorithm capable of stop-and-go. """
  __metaclass__ = ABCMeta
  SHELVECAR = 'SHELVECAR'
  """ Name of the file where the history of a GA run is saved. """
  CALCDIR  = 'calcs'
  """ Name of subdirectory to calculations. """

  @abstractmethod
  def random_individual(self, **kwargs):
    """ Returns new random individual. """
    pass
  @abstractmethod
  def jobinator(self, job, indiv):
    """ Fills jobdictionary item for computing individual.
    
        :param job: 
          :py:class:`JobDict <lada.jobs.JobDict>` instance which should be
          filled with the right parameters for computing the intermediate
          fitness of an individual.
        :param indiv:
          The individual for which to perform a calculation.
    """
    pass
  @abstractmethod
  def objective(self, extract, indiv):
    """ Returns fitness of an individual.
    
        :param extract: 
          An extraction object which comes as a result of the computation
          created by :py:meth:`jobinator <StepGA.jobinator>`. It should be used
          to actually evaluate the fitness of the individual.
        :param indiv:
          Individual evaluate.
    """
    pass

  def update_fitness(self, population=None):
    """ Updates fitness of whole population. 
    
        This is used to update the fitness of the whole population when it
        changes depending on history, e.g. for X-GSGO.
    """
    return 

  def __init__(self, directory=None):
    """" Initializes StepGA. """
    from os.path import exists
    from shelve import open as shelve_open
    from .standard import Mating
    from ..opt import RelativeDirectory
    super(StepGA, self).__init__()

    self._directory = RelativeDirectory(directory)

    self.population = []
    """ List of individuals in the population. """
    self.offspring = []
    """ List of uncomputed individuals. """
    self.generation = 0
    """ Current generation.

        The counter is incremented after each successful call to :py:meth:`next`.
    """
    self.matingops = Mating(sequential=False)
    """ Mating operators. """

    # creates a SHELVECAR file.
    if exists(self.shelvepath):
      raise IOError('{0} already exists.\n'\
                    'Delete it to restart a GA from scratch.\n'\
                    'Use lada.ga.stepped.load() to restart from previous run.'\
                    .format(self.shelvepath))
    try: 
      shelvecar = shelve_open(self.shelvepath)
      shelvecar['individuals'] = {}
      shelvecar['functionals'] = [self]
      shelvecar['added']       = []
      shelvecar['removed']     = []
      shelvecar['new']         = []
    finally: shelvecar.close()

  def selection(self, dummy, size=2):
    """ Selects individual from the population through deterministic tournament. """
    import random
    list_ = range(len(self.population))
    random.shuffle(list_)
    list_ = list_[:size]
    result = list_[0]
    for b in list_[1:]:
      if self.population[b].fitness < self.population[result].fitness: result = b;
    return result

  @property
  def directory(self):
    """ Directory where output should be found. """
    return self._directory.path
  @property
  def shelvepath(self):
    """ Directory where output should be found. """
    from os.path import join
    return join(self._directory.path, self.SHELVECAR)

  @property
  def popsize(self): 
    """ Size of current population. """
    return len(self.population)
  @property
  def nboffspring(self): 
    """ Size of current population. """
    return len(self.offspring)

  def _add_indiv(self, individual, mechanics='mating'):
    """ Adds new individual if it is not taboo. 
    
        Does not update the SHELVECAR file.
    """
    from uuid import uuid4
    if self.is_taboo(individual): return False

    individual.uuid = uuid4()
    individual.conception = self.generation
    individual.mechanics  = mechanics
    self.offspring.append(individual)
    print '{0}: {1}'.format(mechanics, individual)
    return True

  def add_indiv(self, individual, mechanics='manual'):
    """ Adds new individual if it is not taboo. 

        Updates SHELVECAR file. 
    """ 
    if not self._add_indiv(individual, mechanics): return False
    self._update_func()
    return True

  def is_taboo(self, individual):
    """ Returns True if individual is taboo. 


        Generic function checks that the individual is not in offspring or in
        population.
    """
    return individual in self.offspring or individual in self.population

  def mating(self, n=1, itermax=10):
    """ Returns new individual from population. """
    # sanity check
    if len(self.population) == 0:
      print "No population yet..."
      print "Do ga.next() then call this function again"
      return


    dummy = list(self.offspring)
    try:
      # creates new individuals
      for i in xrange(n):
        success = False
        for j in xrange(itermax):
          if self._add_indiv( self.matingops(self) ): 
            success = True
            break
        if not success: raise RuntimeError('Could not create new individual.')
    except: # revert changes on failure.
      self.offspring = dummy
      raise

    self.current()
    self._update_func()


  def new_random_individuals(self, n=1, itermax=10, **kwargs):
    """ Creates a new individual and adds it to the offspring population. """
    # creates new individuals
    dummy = list(self.offspring)
    try:
      for i in xrange(n):
        success = False
        for j in xrange(itermax):
          if self._add_indiv( self.random_individual(**kwargs), mechanics='random' ): 
            success = True
            break
        if not success: raise RuntimeError('Could not create new individual.')
    except: # revert changes on failure.
      self.offspring = dummy
      raise

    self.current()
    self._update_func()

  @property
  def _jobdict(self):
    """ Creates jobdict from offspring. """
    from ..jobs import JobDict
    rootjobdict = JobDict()
    for indiv in self.offspring:
      job = rootjobdict / str(indiv.conception) / str(indiv.uuid)
      job.uuid = indiv.uuid
      job._indiv = indiv
      self.jobinator(job, indiv)
    return rootjobdict

  @property
  def _jobdict_path(self):
    """ Returns current jobdict path. """
    from os.path import join
    return join(self.directory, join(self.CALCDIR, str(self.generation))) + '.dict'

  def current(self, jobdict=None):
    """ Sets ipython jobdict to current offspring. """
    from IPython.ipapi import get as get_ipy
    from ..jobs import MassExtract as Collect, JobParams, save
    ipy = get_ipy()
    ipy.user_ns['current_jobdict'] = self._jobdict if jobdict is None else jobdict
    ipy.user_ns['current_jobdict'].isgadict = True
    if 'collect' not in ipy.user_ns: ipy.user_ns['collect'] = Collect(dynamic=True)
    else: ipy.user_ns['collect'].uncache()
    if 'jobparams' not in ipy.user_ns: ipy.user_ns['jobparams'] = JobParams()
    ipy.user_ns['current_jobdict_path'] = self._jobdict_path
    try: save(ipy.user_ns['current_jobdict'], self._jobdict_path, timeout=10, overwrite=True)
    except RuntimeError: 
      print "Could not acquire lock on jobdictionary."
      print "The jobdictionary was not saved to disk."
      return 
    print "Setting current job-dictionary to GA offspring."

  def next(self, remove=0):
    """ Goes to next GA iteration. """
    from operator import attrgetter, itemgetter
    from shelve import open as shelve_open
    from itertools import chain
    from IPython.ipapi import get as get_ipy
    # sanity check
    if remove > 0 and remove > len(self.population):
      print "Requested to remove {0} individuals from the population,\n"\
            "but the population size is {1}.".format(remove, len(self.population)) 
      return
    # retrieves current jobdictionary if it exists, otherwise creates new one.
    ipy = get_ipy()
    jobdict = ipy.user_ns.get('current_jobdict', None)
    is_current = jobdict is not None
    if is_current == False or 'isgadict' not in jobdict.root.__dict__:
      jobdict = self._jobdict
    else: jobdict = jobdict.root
    if len(jobdict.keys()) == 0:
      print "No offspring, nothing to do."
      return

    # deletes from the offspring those individuals that are not present in the current jobdict.
    if is_current: 
      # check that all jobs in current jobdict are in the population.
      if len(   set([job.uuid for job in jobdict.itervalues()])\
              - set([indiv.uuid for indiv in self.offspring])  ) != 0:
        print "Found job in current jobdictionary which is not part of the offspring. Aborting."
        return

      uuids = set([indiv.uuid for indiv in self.offspring]) - set([job.uuid for job in jobdict.itervalues()])
      for uuid in uuids: 
        for i, indiv in enumerate(self.offspring):
          if indiv.uuid == uuid: break
        self.offspring.pop(i)
   
    # finds out which offspring have been successfully calculated.
    successfuls = {}
    offmap = {}
    for i, indiv in enumerate(self.offspring): offmap[indiv.uuid] = i
    for name, extractor in ipy.user_ns['collect'].iteritems(): 
      if extractor.success == False: continue

      # try and compute that individual's fitness
      uuid = ipy.user_ns['current_jobdict'][name].uuid
      indiv = self.offspring[offmap[uuid]]
      try: indiv.fitness = self.objective(extractor, indiv)
      except:
        print "Error when computing fitness of {0}.".format(uuid)
        raise
      else: successfuls[name] = offmap[uuid]
    if len(successfuls) == 0 and len(self.offspring) == 0: 
      print 'No successful jobs and no new jobs found. Aborting.'
      return

    # update fitness of whole population, in case it depends upon whole population.
    self.update_fitness(chain((self.offspring[i] for i in successfuls.itervalues()), self.population))
  

    # Now reduces population by requested amount if any.
    if len(successfuls) != 0:
      offspring = list(self.offspring)
      if len(self.population) == 0: population, removed = [], [] # at start
      else: # not at start.
        # Minimizing: remove individuals with lowest fitness.
        population = sorted(self.population, key=attrgetter('fitness'))
        # now do removal.
        if remove > 0:
          removed = [u.uuid for u in population[len(population)-remove:]]
          population = population[:len(population)-remove]
        else:
          removed = [u.uuid for u in population[len(population)-len(successfuls):]]
          population = population[:len(population) - len(successfuls)]
      # append/remove successfull individuals.
      for uuid, index in sorted(successfuls.iteritems(), key=itemgetter(1))[::-1]:
        population.append(offspring.pop(index))
        population[-1].birth = self.generation
  
      # check with user. 
      a = 'o'
      while a not in ['n', 'y']:
        a = raw_input( "Removing {0} individual from population.\n"\
                       "Adding {1} offspring to population.\n"\
                       "New population size: {2}\n"\
                       "Remaining offspring (unfinished calculations): {3}\n"
                       "Is this OK? [y/n] "\
                       .format(len(removed), len(successfuls), len(population), len(offspring)))
      if a == 'n': print "Aborted."; return


      # updates population.
      self.population = population
      self.offspring  = offspring
      self.generation += 1

    # update jobdict, disable older individuals, save it.
    # or remove anything to do with jobdictionary if nothing left to compute
    jobdict = self._jobdict
    if len(jobdict.root.keys()) > 0:
      offmap = {}
      for i, indiv in enumerate(self.offspring): offmap[indiv.uuid] = i
      for job in jobdict.itervalues():
        if self.offspring[ offmap[job.uuid] ].conception == self.generation: job.untag()
        else: job.tag()
      self.current(jobdict)
      print "Updated current job-dictionary."
    else:
      ipy.user_ns.pop('current_jobdict', None)
      ipy.user_ns.pop('current_jobdict_path', None)
      ipy.user_ns.pop('collect', None)
      print "No uncomputed jobs left in offspring."


    # save new step into shelf.
    # or simply update SHELVECAR if no successful jobs.
    if len(successfuls) != 0:
      try:
        shelvecar = shelve_open(self.shelvepath, writeback=True)
        individuals = shelvecar['individuals']
        for individual in self.population:  individuals[str(individual.uuid)] = individual
        shelvecar['functionals'].append(self)
        shelvecar['removed'].append(removed)
        shelvecar['added'].append(successfuls.keys())
        shelvecar['new'].append( [ u.uuid for u in self.offspring\
                                   if u.conception == self.generation-1] )
      finally:
        shelvecar.close()
        print "Updated SHELVECAR file."
    else:
      print "No new completed calculation found."
      try:
        shelvecar = shelve_open(self.shelvepath, writeback=True)
        shelvecar['functionals'][-1] = self
      finally:
        shelvecar.close()
        print "Updated SHELVECAR file."

    # tell user to launch on his own.
    if len(jobdict.root.keys()) > 0:
      print "\nUse %launch to run calculation for new individuals."
      print "Note: to also relaunch calculations for older individuals, first turn "\
            "them on with jobparams.onoff = 'on'."

  def _update_func(self):
    """ Updates functional in SHELVECAR file. """
    from shelve import open as shelve_open
    try:
      shelvecar = shelve_open(self.shelvepath, writeback=True)
      shelvecar['functionals'][-1] = self
    finally: shelvecar.close()


def load(): 
  """ Loads StepGA from current directory. """
  from IPython.ipapi import get
  ip = get()
  try: ip.user_ns['extract'] = Extract()
  except:
    print "No SHELVECAR file."
    print "First create a StepGA functional, then call this function again."
    return
  ip.user_ns['ga'] = ip.user_ns['extract'].functionals[-1]
  if len(ip.user_ns['ga'].offspring) > 0: ip.user_ns['ga'].current()

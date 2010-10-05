from ....opt.decorators import make_cached, broadcast_result
import re

__all__  = ['Extract']

class Extract(object): 
  """ Extracts from ga ouput. """
  ordinals = ['first', 'second', 'third', 'fourth', 'fith', 'sixth',
              'seventh', 'eight', 'nineth', 'eleventh', 'twelfth',
              'thirteenth', 'fourteenth', 'fifteenth', 'sixteenth',
              'seventeenth', 'eighteenth', 'nineteenth' ]
  """ Names of each historical age in the GA. """
  OUTCAR = "out"
  """ Output filename. """
  ERRCAR = "err"
  """ Error filename. """
  FUNCCAR = "GA_FUNCCAR"
  """ Functional filename """

  INDIV_PATTERN = re.compile("^\s*(array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*"\
                             "\s*[0,1]\s*\]\))\s*(\S+)\s*(\S+)?\s*$", re.X|re.M)
  OFFSPRING_PATTERN = \
    re.compile("^\s*Offspring:\s*\n"\
               "\s*((array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*\s*[0,1]\s*\]\))"\
               "\s*(\S+)\s*(\S+)?\s*\n)+", re.X|re.M)

  def __init__(self, directory = ".", comm = None):
    """ Initializes Extract object. """
    from lada.opt import RelativeDirectory

    self._directory = RelativeDirectory(path=directory, hook=self.uncache)
    """ GA directory. """
    self.comm = comm
    """ MPI Communicator. """

  def _search_OUTCAR(self, regex, path=None):
    """ Looks for all matches. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    if path == None: 
      assert self.current_age != None, IOError("No data to search. No GA run completed.")
      path = join(join(self.directory, self.current_age), self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    regex  = compile(regex)
    with open(path, "r") as file:
      for line in file: 
        found = regex.search(line)
        if found != None: yield found

  def _find_first_OUTCAR(self, regex, path=None):
    """ Returns first result from a regex. """
    for first in self._search_OUTCAR(regex, path): return first
    return None

  def _rsearch_OUTCAR(self, regex, path=None):
    """ Looks for all matches starting from the end. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    if path == None: 
      assert self.current_age != None, IOError("No data to search. No GA run completed.")
      path = join(join(self.directory, self.current_age), self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    regex  = compile(regex)
    with open(path, "r") as file: lines = file.readlines()
    for line in lines[::-1]:
      found = regex.search(line)
      if found != None: yield found

  def _find_last_OUTCAR(self, regex, path=None):
    """ Returns first result from a regex. """
    for last in self._rsearch_OUTCAR(regex, path): return last
    return None

  @property
  def directory(self):
   """ Root output directory. """
   return self._directory.path
  @directory.setter
  def directory(self, value): self._directory.path = value

  @property
  def success(self):
    """ Was last run successful?
    
       Checks if directories and files exist. 
    """
    from os.path import exists, isdir, join
    
    all_dirs = [ u for u in self.ordinals if exists(join(self.directory, u)) ]
    if len(all_dirs) == 0: return False
    all_dirs = [ u for u in all_dirs if isdir(join(self.directory, u)) ]
    if len(all_dirs) == 0: return len(self.ages) == 0  
    return all_dirs[-1] == self.ages[-1] if len(self.ages) != 0 else False

  @property
  @make_cached
  def functional(self): 
    """ Input functional from last age. """
    from os.path import join, exists
    from pickle import load

    is_mpi = self.comm != None
    is_root = True if not is_mpi else self.comm.rank == 0

    if self.current_age == None: return None
    if is_root:
      current_path = join(join(self.directory, self.current_age), self.FUNCCAR)
      assert exists(current_path), RuntimeError("File {0} does not exist.".format(current_path))
      with open(current_path, "r") as file: result = load(file)
      if is_mpi:
        from boost.mpi import broadcast
        broadcast(self.comm, result, 0)
    elif is_mpi:
      from boost.mpi import broadcast
      result = broadcast(self.comm, None, 0)
    return result

  @property
  def current_age(self):
    """ Current GA age (eg last run). """
    return self.ages[-1] if len(self.ages) else None

  @property
  def next_age(self):
    """ Next GA age (eg next run). """
    if self.current_age == None: return self.ordinals[0] 
    index = self.ordinals.index(self.current_age) + 1
    assert index < len(self.ordinals), RuntimeError("Ran out of ordinals.")
    return self.ordinals[index]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ages(self): 
    """ Returns existing ages. """
    from os.path import exists, join

    results = []
    filenames = [self.OUTCAR, self.FUNCCAR]
    for name in self.ordinals:
      dummy = [join(join(self.directory, name), f) for f in filenames]
      # check for existsnce of all files.
      if not all([exists(p) for p in dummy]): continue
      # checks for the nimber of generations.
      first = self._find_first_OUTCAR( r"^\s*Starting\s+generation\s+(\d+)\s*$",\
                                       join(join(self.directory, name), self.OUTCAR))
      if first == None: continue
      first = int(first.group(1))
      last = self._find_last_OUTCAR( r"^\s*Starting\s+generation\s+(\d+)\s*$",\
                                     join(join(self.directory, name), self.OUTCAR))
      if last == None: continue
      last = int(last.group(1))
      if last - first < 1: continue
      results.append(name)
    return results

  @make_cached
  def solo(self):
    """ Returns extractor with no communicator. """
    from copy import deepcopy
    result = deepcopy(self)
    result.comm = None
    return result

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def start_generation(self):
    """ Generation at start of run. """
    found = self._find_first_OUTCAR("^\s*Starting\s+generation\s+(\d+)\s*$")
    return int(found.group(1)) if found != None else None

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def last_generation(self):
    """ Generation at end of run. """
    found = self._find_last_OUTCAR("^\s*Starting\s+generation\s+(\d+)\s*$")
    return int(found.group(1)) if found != None else None

  @property
  def current_gen(self):
    """ Returns current population. Alias for last_generation. """
    return self.last_generation
  

  @property
  def has_gaps(self):
    """ True if there are gaps in historical ages. """
    if len(self.ages) == 0: return False 
    for a, b in zip(self.ages, self.ordinals):
      if a != b: return True
    return True

  def uncache(self): 
    """ Uncaches results. """
    from lada.opt.decorators  import uncache as opt_uncache
    opt_uncache(self)


  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def offspring(self):
    """ List of bitstring + fitness for each offspring at each generation. """
    from os.path import join
    from numpy import array, abs
    from copy import deepcopy
    
    result = []
    individual = self.functional.Individual()
    for age in self.ages: 
      with open(join(join(self.directory, age), self.OUTCAR), "r") as file:
        text = "".join(file.readlines())

      # loop over Offspring blocks.
      for gen_match in self.OFFSPRING_PATTERN.finditer(text):
        result.append([])
        # loop over inviduals
        for indiv_match in self.INDIV_PATTERN.finditer(gen_match.group(0)):
          # creates individual from grepped stuff.
          individual.genes = eval(indiv_match.group(1))
          individual.fitness = float(indiv_match.group(2))
          # checks for pre-existing copies.
          found_other = False
          if len(result) >= 2: 
            for other in result[-2]: 
              if self.functional.cmp_indiv(individual, other) != 0: continue
              if not self.functional.compare(individual, other): continue
              found_other = True
              break
          result[-1].append( other if found_other else deepcopy(individual) )
      if len(result[-1]) == 0: result.pop(-1)
    return result
               
  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def generations(self):
    """ List of bitstring + fitness for each population at each generation. """
    if len(self.offspring) == 0: return []
    result = [sorted(self.offspring[0], self.functional.cmp_indiv)]
    for generation in self.offspring[1:]:
      new_pop = result[-1][:-len(generation)]
      new_pop.extend(generation)
      result.append( sorted(new_pop, self.functional.cmp_indiv) )
    return result 

  @property
  @broadcast_result(attr=True, which=0)
  def genplots(self):
    """ Generation plotter arrays. """
    from numpy import arange, array, transpose

    if len(self.generations) == 0: return None, None
    y = [ u.fitness for gen in self.generations for u in gen ]
    x = [ i for i, gen in enumerate(self.generations) for u in gen ]

    return x, y

  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def uniques(self):
    """ Minimum set of individuals. """
    results = []
    compare = self.functional.compare
    for pop in self.solo().offspring:
      for indiv in pop:
        i = id(indiv)
        if any( i == id(o) for o in results ): continue
        if any( compare(indiv, o) for o in results ): continue
        results.append(indiv)
    return results

  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def sorted_unique(self):
    """ Minimum set of individuals. """
    from operator import attrgetter
    return sorted(self.solo().uniques, key=attrgetter('fitness'))

  def solo(self):
    """ Returns a serial extractor (as opposed to MPI). """
    from copy import deepcopy
    result = deepcopy(self)
    result.comm = None
    return result

  def __getstate__(self):
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)

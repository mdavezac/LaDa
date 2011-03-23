from ....opt.decorators import make_cached, broadcast_result
from ....opt import AbstractExtractBase, OutcarSearchMixin
import re

__all__  = ['Extract']

class Extract(AbstractExtractBase, OutcarSearchMixin): 
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
  OFFCAR = "GA_OFFCAR"
  """ Offspring file. """

  INDIV_PATTERN = re.compile(r"^\s*(array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*"\
                             r"\s*[0,1]\s*\]\))\s*(\S+)\s*(\S+)?\s*$", re.X|re.M)
  OFFSPRING_PATTERN = \
    re.compile(r"^\s*Offspring:\s*\n"\
               r"(\s*(array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*\s*[0,1]\s*\]\))"\
               r"\s*(\S+)\s*(\S+)?\s*\n)+", re.X|re.M)

  def __init__(self, directory = None, comm = None):
    """ Initializes Extract object. """
    AbstractExtractBase.__init__(self, directory, comm)
    OutcarSearchMixin.__init__(self)

  def __outcar__(self):
    """ Returns path to OUTCAR file.

        :raise IOError: if the OUTCAR file does not exist. 
    """
    from os.path import exists, join
    if self.current_age == None: raise IOError("No data to search. No completed GA run found.")
    path = join(join(self.directory, self.current_age), self.OUTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

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

    if self.current_age == None: return None
    if self.comm.is_root:
      current_path = join(join(self.directory, self.current_age), self.FUNCCAR)
      assert exists(current_path), RuntimeError("File {0} does not exist.".format(current_path))
      with open(current_path, "rb") as file: result = load(file)
      self.comm.broadcast(result)
    else: result = self.comm.broadcast()
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
    from re import compile, M as mult
    from os.path import exists, join

    results = []
    filenames = [self.OUTCAR, self.FUNCCAR]
    for name in self.ordinals:
      dummy = [join(join(self.directory, name), f) for f in filenames]
      # check for existsnce of all files.
      if not all([exists(p) for p in dummy]): continue
      # checks for the nimber of generations.
      regex = compile(r"^\s*Starting\s+generation\s+(\d+)\s*$", mult)
      with open(join(join(self.directory, name), self.OUTCAR), 'r') as file:
        string = file.read()
      first, last = None, 0
      for found in re.finditer(regex, string):
        if first == None: first = int(found.group(1)); continue
        last = int(found.group(1))
      if first == None: continue
      if last - first < 1: continue
      results.append(name)
    return results

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

  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def offspring(self):
    """ List of bitstring + fitness for each offspring at each generation. """
    from copy import deepcopy
    

    def loop(_age):
      from os.path import join, exists
      from numpy import array
      from quantities import *

      path = join(join(self.directory, _age), self.OFFCAR)
      if exists(path):  # loads from OFFCAR
        with open(path, "r") as file: 
          for pop in load(file):
            for indiv in pop: yield indiv
      else:
        individual = self.functional.Individual()
        with open(join(join(self.directory, _age), self.OUTCAR), "r") as file:
          text = "".join(file.readlines())
        
        # loop over Offspring blocks.
        for gen_match in self.OFFSPRING_PATTERN.finditer(text):
          result.append([])
          # loop over inviduals
          for indiv_match in self.INDIV_PATTERN.finditer(gen_match.group(0)):
            individual.genes = eval(indiv_match.group(1))
            individual.fitness = float(indiv_match.group(2))
            yield individual

    result = []
    for age in self.ages: 
      for individual in loop(age):
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

  def uniquify_lists(self, *args): 
    """ Make list of unique individuals from different lists. """
    from itertools import chain
    results = []
    compare = self.functional.compare
    for indiv in chain(*args):
      i = id(indiv)
      if any( i == id(o) for o in results ): continue
      if not any( compare(indiv, o) for o in results ): results.append(indiv)
    return results

  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def uniques(self):
    """ Minimum set of individuals. """
    return self.uniquify_lists(*self.solo().offspring)

  @property
  @make_cached 
  @broadcast_result(attr=True, which=0)
  def sorted_unique(self):
    """ Minimum set of individuals. """
    from operator import attrgetter
    return sorted(self.solo().uniques, key=attrgetter('fitness'))

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        :kwarg errors: Include stderr files.
        :type errors: bool
        :kwarg offcar: Include file with offspring for each generation.
        :type offcar: bool
    """
    from os.path import join, exists, isdir
    files = [self.OUTCAR, self.FUNCCAR]
    if kwargs.get("errors", False): files.append(self.ERRCAR)
    if kwargs.get("offcar", False): files.append(self.OFFCAR)

    for age in self.ages:
      directory = join(self.directory, age)
      if not exists(directory): continue
      if not isdir(directory): continue
      for file in files:
        path = join(directory, file)
        if exists(path): yield path
    


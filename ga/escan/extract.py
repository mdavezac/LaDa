from lada.opt.decorators import make_cached 
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
  INCAR = "input.py"
  """ Input filename. """
  FUNCCAR = "oscillator.py"
  """ Functional filename. """

  INDIV_PATTERN = re.compile("^\s*(array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*"\
                             "\s*[0,1]\s*\]\))\s*(\S+)\s*$", re.X|re.M)
  OFFSPRING_PATTERN = \
    re.compile("^\s*Offspring:\s*\n"\
               "\s*((array\(\[(?:\s*[0,1]\,|(?:\,\s*$^)\s*)*\s*[0,1]\s*\]\))"\
               "\s*(\S+)\s*\n)+", re.X|re.M)

  def __init__(self, directory = "."):
    """ Initializes Extract object. """
    from lada.opt import RelativeDirectory

    root = RelativeDirectory(path=directory, hook=lambda x:self.uncache())
    """ GA directory. """


  @property
  @make_cached
  def input(self): 
    """ Input parameters. """
    from os.path import join, exists
    from lada.ga.escan.elemental import Converter
    from lada.ga.escan.elemental.evaluator import Dipole 
    from lada.opt   import read_input
    from lada.vff   import Vff
    from lada.escan import Escan, soH

    assert exists(join(self.root, self.OUTCAR)), \
           IOError("Could not find file {0}.".format(join(self.root, self.OUTCAR)))
    return read_input\
           (\
             join(self.root.path, self.OUTCAR),\
             { "Converter": Converter, "DipoleEvaluator": Dipole,\
               "Vff": Vff, "Escan": Escan, "soH": soH },\
             paths = ("saveto_path", "restartfrom_path", "workdir", "outdir")
           )

  @property
  @make_cached
  def functional(self): 
    """ Input parameters. """
    from os.path import join, exists
    from lada.opt import read_input

    assert exists(join(self.root.path, self.FUNCCAR)), \
           IOError("Could not find file {0}.".format(join(self.root.path, self.FUNCCAR)))
    funccar = read_input(join(self.root.path, self.FUNCCAR) )
    return funccar.Darwin()

  @property
  @make_cached
  def Functional(self): 
    """ Returns functional type. """
    from os.path import join, exists
    from lada.opt import read_input

    assert exists(join(self.root.path, self.FUNCCAR)), \
           IOError("Could not find file {0}.".format(join(self.root.path, self.FUNCCAR)))
    funccar = read_input(join(self.root.path, self.FUNCCAR) )
    return funccar.Darwin
   
  @property
  @make_cached
  def ages(self): 
    """ Returns existing ages. """
    from os.path import exists, join

    results = []
    for name in self.ordinals:
      filepath = join(join(self.root.path, name), self.OUTCAR)
      filepath = join(join(self.root.path, name), self.FUNCCAR)
      if exists(filepath): results.append(name)
    return results

  @property
  def has_gaps(self):
    """ True if there are gaps in historical ages. """
    if len(self.ages) == 0: return False 
    for a, b in zip(self.ages, self.ordinals):
      if a != b: return True
    return True

  @property
  def next_dirname(self):
    """ Returns next directory name. """
    return self.ordinals[self.ordinals.index(self.ages[-1]) + 1]

  def uncache(self): 
    """ Uncaches results. """
    from lada.opt.decorators  import uncache as opt_uncache
    opt_uncache(self)

  def move_to_new_age(self): 
    """ Moves current results to a new age directory. """
    from os import makedirs
    from os.path import join, exists
    from shutil import copy, move

    # uncache results, make sure we don't overwrite anything.
    self.uncache()

    if not exists(join(self.root.path, "stdout")): return
    newdir = join(self.root.path, self.next_dirname)
    if not exists(newdir): makedirs(newdir)
    for file in ["stdout", "stderr"]:
      move(join(self.root.path, file), join(newdir, self.OUTCAR))
    for file in [self.INCAR, self.FUNCCAR]: 
      copy(join(self.root.path, file), newdir)


  @property
  @make_cached 
  def offspring(self):
    """ List of bitstring + fitness for each offspring at each generation. """
    from os.path import join
    from numpy import array, abs
    from copy import deepcopy
    
    result = []
    individual = self.functional.Individual()
    for age in self.ages: 
      with open(join(join(self.root.path, age), self.OUTCAR), "r") as file:
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
              if any( abs(individual.genes - other.genes[0]) > 1e-12 ): continue
              found_other = True
              break
          result[-1].append( other if found_other else deepcopy(individual) )
      if len(result[-1]) == 0: result.pop(-1)
    return result
               
  @property
  @make_cached 
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
  def genplots(self):
    """ Generation plotter arrays. """
    from numpy import arange, array, transpose

    if len(self.generations) == 0: return None, None
    y = [ u.fitness for gen in self.generations for u in gen ]
    x = [ i for i, gen in enumerate(self.generations) for u in gen ]

    return x, y

  def __getstate__(self):
    d = self.__dict__.copy()
    d.pop("comm", None)
    return d
  def __setstate__(self, arg):
    self.__dict__.update(arg)

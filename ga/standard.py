#
#  Version: $Id$
#
""" Holds standard genetic algorithm operations. """
def no_taboo( self, _indiv ):
  """ Does nothing. """
  return False
def taboo( self, _indiv ):
  """ taboo makes sure that no two individuals in the population and the offspring are the same. """
  return _indiv in self.population or _indiv in self.offspring

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
  """ Compares two individuals. """
  from math import fabs
  if fabs(a.fitness - b.fitness) <  tolerance: return 0
  elif a.fitness < b.fitness: return -1
  return 1

def average_fitness(self):
  """ Prints out average fitness. """
  result = 0e0
  for indiv in self.population:
    result += indiv.fitness
  print "  Average Fitness: ", result / float(len(self.population))
  return True

def best(self):
  best = None
  for indiv in self.population:
    if best == None or  self.cmp_indiv( best, indiv ) == 1: 
      best = indiv
  print "  Best Individual: ", best.genes, best.fitness
  return True

def print_population(self):
  print "  Population: "
  for indiv in self.population:
    print "    ", indiv.genes, indiv.fitness
  return True

def print_offspring(self):
  print "  Offspring: "
  for indiv in self.population:
    if indiv.birth == self.current_gen - 1: 
      print "    ", indiv.genes, indiv.fitness
  return True

def check_generation( self ):
  """ returns false if maximum number of generations was passed. """
  if self.max_gen < 0: return True
  return self.current_gen < self.max_gen
  
def add_population_evaluation(self, evaluation):
  """ Standard population evaluation. 
      Evaluates individual only if fitness attribute does not exist. 
      Fitness is the return of evaluation subroutine given on input.
      evaluation subroutine should take an individual at its only argument.
  """
  def popeval(self):
    for indiv in self.population:
      if not hasattr(indiv, "fitness" ): 
        indiv.fitness = self.indiv_evaluation(indiv)
    for indiv in self.offspring:
      if not hasattr(indiv, "fitness" ): 
        indiv.fitness = self.indiv_evaluation(indiv)
  self.indiv_evaluation = evaluation
  self.evaluation = popeval
  return self

def add_checkpoint(self, _chk):
  """ Adds a checkpoint """
  try: self.checkpoints.append( _chk ) 
  except AttributeError: self.checkpoints = [_chk]
  return self;

def fill_attributes(self):
  """ Checks self for correct attributes.
      Fills in where possible:
        _ "taboo" defaults to standard.taboo
        _ "selection" defaults to standard.selection
        _ "population" defaults to empty list []
        _ "popsize" defaults to len(self.population)
        _ "offspring" defaults to empty list []
        _ "rate" defaults to len(self.offspring)/self.popsize
        _ standard.check_generation is ALWAYS added to the checkpoints
        _ "current_gen" defaults to 0 
        _ "max_gen" defaults to 100, or current_gen+100 if max_gen > current_gen.
        _ "cmp_indiv" defaults to standard.cmp_indiv
  """
  import darwin 
  # must have an evaluation function.
  assert hasattr(self, "evaluation"), "No evaluation function!" 

  # Checks that self has an object Individual
  if not hasattr(self, "Individual"):
    from bitstring import Individual
    self.Individual = Individual

  # Checks whether self has a taboo object.
  if not hasattr(self, "taboo"): self.taboo = taboo

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
  self = add_checkpoint(self, check_generation)

  # checks current generation.
  if not hasattr(self, "current_gen"): self.current_gen = 0
  elif not hasattr(self, "max_gen"): self.max_gen = 100
  elif self.max_gen < self.current_gen: self.max_gen += 100

  if not hasattr(self, "max_gen"): self.max_gen = self.current_gen + 100

  if not hasattr(self, "cmp_indiv"): self.cmp_indiv = cmp_indiv
  return self

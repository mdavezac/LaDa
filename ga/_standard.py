#
#  Version: $Id$
#
def no_taboo( self, _indiv ): return False
def taboo( self, _indiv ): return _indiv in self

def tournament( self, population, _size = 2 ):
  """ deterministic tournament """
  import random
  list_ = random.shuffle( range(len(population)) )[:_size]
  result = 0
  for b in list_:
    if population[b].fitness < population[result].fitness: result = b;
  return result

def check_generation( darwin ):
  """ returns false if maximum number of generations was passed. """
  if darwin.max_gen < 0: return True
  return darwin.current_gen < darwin.max_gen

def add_population_evaluation(darwin, evaluation):
  """ Standard population evaluation. """
  def popeval(self):
    for indiv in self.population:
      if not hasattr( indiv.fitness ): 
        indiv.fitness = self.indiv_evaluation(indiv)
  darwin.indiv_evaluation = evaluation
  darwin.evaluate = popeval

  
def fill_darwin(darwin):
  """ Checks darwin for correct attributes.
      Fills in where possible. 
  """
  # must have an evaluation function.
  assert hasattr(darwin, "evaluation"), "No evaluation function!" 

  # Checks that darwin has an object Individual
  if not hasattr(darwin, "Individual"):
    from bitstring import Individual
    add_individual_type(darwin, Individual)

  # Checks whether darwin has a taboo object.
  if not hasattr(darwin, "taboo"): add_taboo(darwin, taboo)

  # Checks whether darwin has a selection object.
  if not hasattr(darwin, "selection"): add_selection(darwin, tournament)


  # checks whether there is a population.
  if not hasattr(darwin, "population"): darwin.population = []
  
  # checks whether there is a popsize.
  if not hasattr(darwin, "popsize"): darwin.popsize = len(darwin.population)

  # makes sure we have something to do.
  assert darwin.popsize != 0, "population or popsize attributes required on input."

  # checks whether there is a population.
  if not hasattr(darwin, "offspring"): darwin.offspring = []
  
  # checks whether there is a popsize.
  if not hasattr(darwin, "rate"):
    darwin.rate = float(len(darwin.offspring)) / float(darwin.popsize))

  # makes sure we have something to do.
  assert darwin.rate > float(0), "offspring or rate attributes required on input."
 
  # checks whether there is a checkpoint.
  add_checkpoint(darwin, check_generation)

  # checks current generation.
  if not hasattr(darwin, "current_gen"): darwin.current_gen = 0
  elif not hasattr(darwin, "max_gen"): darwin.max_gen = 100
  elif darwin.max_gen < darwin.current_gen: darwin.max_gen += 100

  if not hasattr(darwin, "max_gen"): darwin.max_gen = darwin.current_gen + 100

  return darwin

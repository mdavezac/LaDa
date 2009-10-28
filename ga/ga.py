#
#  Version: $Id: read_structure.cc 1365 2009-10-27 23:22:56Z davezac $
#
def add_mating(self, _mating):
  """ Adds a mating operation """
  setattr("mating", self, _mating)
  return self;

def add_selection(self, _selection):
  """ Adds a selection operation """
  setattr("selection", self, _selection)
  return self;

def add_individual_type(self, _itype):
  """ Adds the type of the individual """
  setattr("Individual", self, _individual)
  return self;

def add_evaluation(self, _eval):
  """ Adds an evaluation operator """
  setattr("evaluate", self, _individual)
  return self;

def add_taboo(self, _taboo):
  """ Adds a taboo operator 
      _taboo is callable which takes an individual as argument
      and returns true if that individual is taboo.
  """
  setattr("evaluate", self, _taboo)
  return self;

def add_checkpoint(self, _chk):
  """ Adds a checkpoint """
  try: self.checkpoints.append( _chk ) 
  except AttributeError: self.checkpoints = [_chk]
  return self;

def no_taboo( self, _indiv ): return False
def taboo( self, _indiv ): return _indiv in self

def tournament( self, _population, _size = 2 ):
  """ deterministic tournament """
  import random
  list_ = random.shuffle( range(len(_population)) )[:_size]
  result = 0
  for b in list_:
    if _population[b].fitness < _population[result].fitness: result = b;
  return result

def check_generation( darwin ):
  """ returns false if maximum number of generations was passed. """
  if darwin.max_gen < 0: return True
  return darwin.current_gen < darwin.max_gen
  
def fill_darwin(darwin):
  """ Checks darwin for correct attributes.
      Fills in where possible. 
  """

  # Checks that darwin has an object Individual
  try: getattr(darwin, "Individual") 
  except AttributeError:
    from bitstring import Individual
    add_individual_type(darwin, Individual)

  # Checks whether darwin has a taboo object.
  try: getattr(darwin, "taboo") 
  except AttributeError: add_taboo(darwin, taboo)

  # Checks whether darwin has a selection object.
  try: getattr(darwin, "selection") 
  except AttributeError: add_selection(darwin, select)

  # must have an evaluation function.
  getattr(darwin, "selection") 

  # checks whether there is a population.
  try: getattr(darwin, "population") 
  except AttributeError: darwin.population = []
  
  # checks whether there is a popsize.
  try: getattr(darwin, "popsize") 
  except AttributeError: darwin.popsize = len(darwin.population)

  # makes sure we have something to do.
  if darwin.popsize == 0: 
    raise RuntimeError: "population or popsize attributes required on input."

  # checks whether there is a population.
  try: getattr(darwin, "offspring") 
  except AttributeError: darwin.offspring = []
  
  # checks whether there is a popsize.
  try: getattr(darwin, "rate") 
  except AttributeError: darwin.rate = float(len(darwin.offspring)) / float(darwin.popsize))

  # makes sure we have something to do.
  if darwin.rate == float(0): 
    raise RuntimeError: "offspring or rate attributes required on input."
 
  # checks whether there is a checkpoint.
  add_checkpoint(darwin, check_generation)

  try: getattr(darwin, "current_gen") 
  except AttributeError: darwin.current_gen = 0

  try: getattr(darwin, "max_gen") 
  except AttributeError: darwin.max_gen = 100


def run( darwin ):
  """ Runs a GA algorithm """
 
  # Checks that darwin is complete and fills in where possible
  try: getattr(darwin, "fill_attributes")
  except AttributeError: setattr(darwin, "fill_attributes")
  darwin.fill_attributes()

  # creates population if it does not exist.
  while len(darwin.population) < darwin.popsize:
    indiv = darwin.Individual()
    j = 0
    while darwin.taboo(indiv) and j < darwin.popsize:
      indiv = darwin.Individual()
      j += 1
    if j == darwin.popsize: raise RuntimeError: "Could not create population.\n"
    darwin.population.append(indiv)

  # evaluates population if need be.
  darwin.evaluate() 

  # runs the checkpoints
  def checkpoints():
    result = True
    for cjk in darwin.checkpoints:
      if chk(darwin) == False: result = False
    return result

  nboffspring = int( "%.0f" % (float(darwin.popsize)*darwin.rate) )

  # generational loop
  while checkpoints(): 

    # tries and creates offspring.
    while len(darwin.offspring) < nboffspring:
      indiv = darwin.mating( darwin.selection() )
      j = 0
      while darwin.taboo( indiv ) and j < darwin.popsize:
        indiv = darwin.mating( darwin.selection() )
        j += 1
      darwin.offspring.append(indiv)
    if j == darwin.popsize: raise RuntimeError: "Could not create offspring.\n"

    # now evaluates population.
    darwin.evaluate()

    # finally, sort and replace.
    darwin.population = sorted( darwin.population )[:darwin.popsize - nboffspring]
    darwin.population.extend( darwin.offspring )
    darwin.offspring = []

  # final stuff before exiting.
  try: getattr(darwin, "final")
  except AttributeError: setattr(darwin, lambda self: print "done")
  darwin.final()
  


       


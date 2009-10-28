#
#  Version: $Id$
#
import standard
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

def run( darwin ):
  """ Runs a GA algorithm """
 
  # Checks that darwin is complete and fills in where possible
  if hasattr(darwin, "fill_attributes"): darwin = darwin.fill_attributes()
  else: darwin = standard.fill_attributes(darwin)

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
        indiv = darwin.mating( darwin.selection, darwin.population )
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
  


       


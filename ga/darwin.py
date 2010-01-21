""" Holds a single function for running a ga algorithm """
def run(self):
  """ Performs a Genetic Algorithm search """
  import standard
 
  # runs the checkpoints
  def checkpoints():
    result = True
    for chk in self.checkpoints:
      if chk(self) == False: result = False
    return result


  # Checks that self is complete and fills in where possible
  if hasattr(self, "fill_attributes"): self = self.fill_attributes()
  else: self = standard.fill_attributes(self)
  assert hasattr(self, "Individual"), "No Individual type.\n"
  assert hasattr(self, "taboo"), "No taboo operation.\n"
  assert hasattr(self, "population"), "No population attribute.\n"
  assert hasattr(self, "popsize"), "No popsize attribute.\n"
  assert hasattr(self, "current_gen"), "No current_gen attribute.\n"
  assert hasattr(self, "evaluation"), "No evaluation operation.\n"
  assert hasattr(self, "checkpoints"), "No checkpoints attribute.\n"
  assert hasattr(self, "rate"), "No rate attribute.\n"
  assert hasattr(self, "mating"), "No mating functor.\n"
  assert hasattr(self, "offspring"), "No offspring attribute.\n"
  assert hasattr(self, "cmp_indiv"), "No cmp_indiv operation.\n"

  # creates population if it does not exist.
  while len(self.population) < self.popsize:
    indiv = self.Individual()
    j = 0
    loop = True
    while loop:
      indiv = self.Individual()
      loop = self.taboo(self, indiv)
      j += 1
      assert j < max(50*self.popsize, 100), "Could not create offspring.\n"
    indiv.birth = self.current_gen
    self.population.append(indiv)

  # evaluates population if need be.
  self.evaluation(self) 

  nboffspring = int( "%.0f" % (float(self.popsize)*self.rate) )

  # generational loop
  while checkpoints(): 
    print "Starting generation ", self.current_gen

    # tries and creates offspring.
    while len(self.offspring) < nboffspring:
      j = 0
      loop = True
      while loop:
        indiv = self.mating(self)
        loop = self.taboo(self, indiv)
        j += 1
        assert j < max(10*self.popsize, 100), "Could not create offspring.\n"
      indiv.birth = self.current_gen
      self.offspring.append(indiv)

    # now evaluates population.
    self.evaluation(self)

    # finally, sort and replace.
    self.population = sorted( self.population, self.cmp_indiv )[:len(self.population)-nboffspring]
    self.population.extend( self.offspring )
    
    self.offspring = []
    self.current_gen += 1 


  # final stuff before exiting.
  if hasattr(self, "final"): self.final()
  else: print "done"
  


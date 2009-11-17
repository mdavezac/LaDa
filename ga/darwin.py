#
#  Version: $Id$
#
def run( self ):
  """ Runs a GA algorithm """
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

  # creates population if it does not exist.
  while len(self.population) < self.popsize:
    indiv = self.Individual()
    j = 0
    loop = True
    while loop:
      indiv = self.Individual()
      loop = self.taboo(self, indiv)
      j += 1
      assert j < self.popsize, "Could not create population. %i\n" % (j)
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
        assert j < self.popsize, "Could not create offspring.\n"
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
  


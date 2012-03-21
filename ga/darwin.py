""" Holds a single function for running a ga algorithm """

def is_synced(comm):
  if not __debug__: return
  result = comm.broadcast("Checking if processes are synced.", 0)
  assert result == "Checking if processes are synced.", \
         RuntimeError("Processes are not synced: %s " % result) 

def _check_population(self, population):
  if not __debug__: return
  new_pop = self.comm.broadcast(population)
  for a, b in zip(new_pop, population):
    assert a == b, RuntimeError("Populations are not equivalent across processes.")
    ahas, bhas = hasattr(a, "fitness"), hasattr(b, "fitness")
    assert (ahas and bhas) or not (ahas or bhas),\
          RuntimeError("Populations do not have equivalently evaluated invidiuals")
    if ahas and bhas:
      assert a.fitness == b.fitness,\
      RuntimeError("Fitness are not equivalent across processes.")
  self.comm.barrier()


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
  if hasattr(self, "fill_attributes"): self.fill_attributes()
  else: self = standard.fill_attributes(self)
  assert hasattr(self, "Individual"), "No Individual type.\n"
  assert hasattr(self, "taboo"), "No taboo operation.\n"
  assert hasattr(self, "population"), "No population attribute.\n"
  assert hasattr(self, "popsize"), "No popsize attribute.\n"
  assert hasattr(self, "current_gen"), "No current_gen attribute.\n"
  assert hasattr(self, "evaluation"), "No evaluation operation.\n"
  assert hasattr(self, "rate"), "No rate attribute.\n"
  assert hasattr(self, "mating"), "No mating functor.\n"
  assert hasattr(self, "offspring"), "No offspring attribute.\n"
  assert hasattr(self, "comparison"), "No comparison operation.\n"

  # creates population if it does not exist.
  if self.comm.is_root: 
    while len(self.population) < self.popsize:
      j = 0
      loop = True
      while loop:
        indiv = self.Individual()
        loop = self.taboo(indiv)
        j += 1
        if j > max(50*self.popsize, 100):
          raise Exception("Could not create initial population.\n")
      indiv.birth = self.current_gen
      self.population.append(indiv)
  self.population = self.comm.broadcast(self.population)
  
  # evaluates population if need be.
  self.evaluation() 
  # now makes sure evaluation did not create twins.
  dummy, self.population = self.population, []
  for indiv in dummy:
    if any(indiv == u for u in self.population): continue
    self.population.append(indiv)

  # Number of offspring per generation.
  nboffspring = max(int(float(self.popsize)*self.rate), 1)

  # generational loop
  while checkpoints(): 
    if self.comm.do_print:
      print "\nStarting generation ", self.current_gen

    # tries and creates offspring.
    if self.comm.rank == 0:
      while len(self.offspring) < nboffspring:
        j = 0
        loop = True
        while loop:
          indiv = self.mating()
          loop = self.taboo(indiv)
          j += 1
          if j > max(10*self.popsize, 100): raise Exception("Could not create offspring.\n")
        indiv.birth = self.current_gen
        self.offspring.append(indiv)
      self.comm.broadcast(self.offspring)
    else: self.offspring = self.comm.broadcast(self.offspring)


    # now evaluates population.
    self.evaluation()

    # finally, sort and replace.
    if self.comm.rank == 0:
      # deal with differences in function sorted between python versions.
      self.population = sorted(self.population, key=self.comparison)
      # Inserts individual one by one ensuring they do not yet exist in the
      # population. This ensures that duplicates are not allowed.
      for indiv in self.offspring:
        if any(indiv == u for u in self.population): continue
        self.population.insert(0, indiv)
        # In case population smaller than expected, use conditional.
        if len(self.population) > self.popsize: self.population.pop(-1)
    self.population = self.comm.broadcast(self.population)
    
    self.offspring = []
    self.current_gen += 1 

  # final stuff before exiting.
  if hasattr(self, "final"): self.final()
  elif self.comm.do_print: print "done"
  


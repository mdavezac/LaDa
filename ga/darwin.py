""" Holds a single function for running a ga algorithm """

def is_synced(comm):
  if not __debug__: return

  from boost.mpi import broadcast
  result = broadcast(comm, "Checking if processes are synced.", 0)
  assert result == "Checking if processes are synced.", \
         RuntimeError("Processes are not synced: %s " % result) 

def _check_population(self, population):
  from boost.mpi import broadcast
  if not __debug__: return
  new_pop = broadcast(self.comm, population, 0)
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
  from boost.mpi import broadcast
  import standard
  import sys
 
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
  assert hasattr(self, "rate"), "No rate attribute.\n"
  assert hasattr(self, "mating"), "No mating functor.\n"
  assert hasattr(self, "offspring"), "No offspring attribute.\n"
  assert hasattr(self, "cmp_indiv"), "No cmp_indiv operation.\n"

  # creates population if it does not exist.
  while len(self.population) < self.popsize:
    j = 0
    loop = True
    while loop:
      indiv = broadcast(self.comm, self.Individual() if self.comm.rank == 0 else None, 0)
      loop = self.taboo(indiv)
      j += 1
      assert j < max(50*self.popsize, 100), "Could not create offspring.\n"
    indiv.birth = self.current_gen
    self.population.append(indiv)

  # evaluates population if need be.
  self.evaluation() 

  # Number of offspring per generation.
  nboffspring = max(int(float(self.popsize)*self.rate), 1)

  # generational loop
  while checkpoints(): 
    if self.comm.do_print:
      print "Starting generation ", self.current_gen

    # tries and creates offspring.
    if self.comm.rank == 0:
      while len(self.offspring) < nboffspring:
        j = 0
        loop = True
        while loop:
          indiv = self.mating()
          loop = self.taboo(indiv)
          j += 1
          assert j < max(10*self.popsize, 100), "Could not create offspring.\n"
        indiv.birth = self.current_gen
        self.offspring.append(indiv)
      broadcast(self.comm, self.offspring, 0)
    else: self.offspring = broadcast(self.comm, self.offspring, 0)


    # now evaluates population.
    self.evaluation()

    # finally, sort and replace.
    if self.comm.rank == 0:
      self.population = sorted( self.population, self.cmp_indiv )[:len(self.population)-nboffspring]
      self.population.extend( self.offspring )
    self.population = broadcast(self.comm, self.population, 0)

    
    self.offspring = []
    self.current_gen += 1 


  # final stuff before exiting.
  if hasattr(self, "final"): self.final()
  elif self.comm.do_print: print "done"
  


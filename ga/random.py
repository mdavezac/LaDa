"""" Holds a single function for running a random search """
def run(self):
  """ Performs a random search """
  import standard
  
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
  assert hasattr(self, "evaluation"), "No evaluation operation.\n"
  assert hasattr(self, "checkpoints"), "No checkpoints attribute.\n"

  while checkpoints():
    indiv = self.Individual()
    j = 0
    loop = True
    while loop:
      indiv = self.Individual()
      loop = self.taboo(self, indiv)
      j += 1
      assert j < max(50*self.popsize, 100), "Could not create offspring.\n"
    self.population = [indiv]
    self.evaluation(self) 

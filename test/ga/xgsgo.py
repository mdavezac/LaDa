from lada.ga import load
def create(input='input.py'):
  """ Load ga into ipython. """
  from os.path import exists
  from pickle import dump
  from lada.vasp import read_input
  from lada.ga.xgsgo import Functional
  from lada.jobfolder import JobFolder

  input = read_input(input)


  root = JobFolder()
  for trial in xrange(input.trials):
    folder = root / str(trial)
    folder.functional = Functional( functional        = input.vasp, 
                                    species           = input.species,
                                    natoms            = input.natoms,
                                    rate              = input.rate,
                                    popsize           = input.popsize,
                                    cns_rate          = input.cns_rate,
                                    mix_atoms_rate    = input.mix_atoms_rate,
                                    mix_poscar_rate   = input.mix_poscar_rate,
                                    mix_mutation_rate = input.mix_mutation_rate )

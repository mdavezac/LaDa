def create(input='input.py'):
  """ Load ga into ipython. """
  from os.path import exists
  from pickle import dump
  from pylada.vasp import read_input
  from pylada.ga.xgsgo.functional import Functional
  from pylada.jobfolder import JobFolder

  input = read_input(input)


  root = JobFolder()
  for trial in xrange(input.trials):
    folder = root / "results" / str(trial)
    folder.functional = Functional( functional        = input.vasp, 
                                    species           = input.species,
                                    natoms            = input.natoms,
                                    rate              = input.rate,
                                    popsize           = input.popsize,
                                    cns_rate          = input.cns_rate,
                                    mix_atoms_rate    = input.mix_atoms_rate,
                                    mix_poscar_rate   = input.mix_poscar_rate,
                                    mutation_rate = input.mutation_rate )
  return root

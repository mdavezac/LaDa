species = ['C', 'N']
""" Species involved in the optimization. """
anions  = ['N']
""" Anions, if any, among the species.

    All anions must be included in the species.
"""
natoms  = (2, 20)
""" Min, Max number of atoms. """
rate    = 0.1
""" Offspring mutation rate.  """
popsize = 50
""" Target population size. """
cns_rate = 0.8
""" Cut and splice crossover rate. """
mutation_rate = 0.1
""" Cut and splice crossover rate. """
mix_atoms_rate = -1
""" mix-atom crossover rate. """
mix_poscar_rate = -1
""" mix-poscar crossover rate. """
trial = 4
""" Number of independent trials. """

first_trial = {'kpoints': '\n0\nAuto\n1\n'}
vasp = Relax(encut=1.0, kpoints='\n0\nAuto\n10', relaxation='volume ionic cellshape')
vasp.species = load('vasp/fere/pseudos', *species)

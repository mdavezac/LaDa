""" Input script to CE-dos GA. """
# vff parameters
from lada.crystal.binary import zinc_blende
lattice = zinc_blende()
for site in lattice.sites: site.type = "Si", "Ge"

directory = "../start"
""" Input directories with calculations. """

clusters = {"J0":True, "J1":True, "B2":12, "B3":4, "B4":2, "B5":2}

# GA parameters.
population_size = 8
""" Size of the GA population """
max_generations = -1
""" Maximum number of generations """
offspring_rate  = 0.25
""" Rate at which offspring are created. """
crossover_rate = 0.75
""" Rate of crossover over other operations. """
trials = range(3)
""" Number of trial GA runs. """
comparison = minimize
""" What we are trying to do with the fitness. """
history = True
""" Whether to keep track of candidates in order never to visit twice. """
histlimit = int(1.5*population_size)
""" Maximum number of individuals to keep in history. """

energies = arange(11, dtype="float64") / 10. * 5 - 2.5
""" Energies over which to perform DOS. """ 
sigma = 0.1 * eV
""" Gaussian smearing for density of states. """

maxiter = -1
""" Maximum of iterations during local search. """
maxeval = 50
""" Maximum of evaluations during local search. """

alwayson = [0, 1]
""" Always leave these figures on. """
alwaysoff = []
""" Always leave these figures off. """

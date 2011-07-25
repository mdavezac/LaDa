""" Input script to CE-dos GA. """
# vff parameters
lattice = crystal.binary.zinc_blende()
for site in lattice.sites: site.type = "Si", "Ge"

directory = "."
""" Input directories with calculations. """

clusters = {"J0":True, "J1":True, "B2":12, "B3":5 "B4":3, "B5":2}

# GA parameters.
population_size = 100
""" Size of the GA population """
max_generations = -1
""" Maximum number of generations """
offspring_rate  = 0.1
""" Rate at which offspring are created. """
crossover_rate = 0.75
""" Rate of crossover over other operations. """
swap_rate = 0.1
""" Rate of swapping mutation over other operations. """
growth_rate = 0.15
""" Rate of growth/shrink mutation over other operations. """
trials = range(3)
""" Number of trial GA runs. """

from numpy import arange
from quantities import eV
energies = arange(11, dtype="float64") / 10. * 5 - 2.5
""" Energies over which to perform DOS. """ 
sigma = 0.1 * eV
""" Gaussian smearing for density of states. """

maxiter = -1
""" Maximum of iterations during local search. """
maxeval = -1
""" Maximum of evaluations during local search. """

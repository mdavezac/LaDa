""" Input script to oscillator strength GA. """
# vff parameters
vff = Vff()
vff.lattice.set_types = ("Si", "Ge"), ("Si","Ge")
vff.add_bond = "Si", "Si", (2.35905923, 45.375785351, -150.363195485, 519.230157133)
vff.add_bond = "Ge", "Ge", (2.44167023, 37.816422666, -124.830189020, 250.179861133)
vff.add_bond = "Si", "Ge", (2.39642606, 41.875535703, -145.994430120, 366.853558523)
vff.add_angle = "Si", "Si", "Si", ("tet", 0, 12.579328566, 44.930324684, 359.422663897)
vff.add_angle = "Ge", "Ge", "Ge", ("tet", 0, 10.383093493, 55.677487481)
vff.add_angle = "Si", "Ge", "Si", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Ge", "Si", "Ge", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Si", "Si", "Ge", ("tet", 0, 12.245971457, 32.561864518, 179.706525419)
vff.add_angle = "Si", "Ge", "Ge", ("tet", 0, 11.147853920, 37.930591733)
vff.minimizer.verbose = True
vff.minimizer.type = "gsl_bfgs2"
vff.minimizer.itermax = 4000
vff.minimizer.tolerance = 1e-5
vff.minimizer.uncertainties = 1e-3
vff.direction = FreezeCell.a1 | FreezeCell.a2

# escan parameters
escan = EMassFunctional()
escan.fft_mesh              = fftmesh.SmallCells()
escan.cutoff                = 8.2
escan.smooth                = 1.0
escan.kinetic_scaling       = 1.0
escan.nbstates              = 10
escan.itermax               = 20
escan.nlines                = 200
escan.tolerance             = 1e-6
escan.rspace_cutoff         = 3.75
escan.do_genpot             = True
escan.do_relax              = True
escan.input_wavefunctions   = None
escan.kpoint                = array([0,0,0])
escan.dnc_mesh              = None
escan.overlap_mesh          = None
escan.potential             = soH
escan.add_potential         = "~/SiGe/pseudos/vq.Si",      "~/SiGe/pseudos/vwr.pso", 0,0,0,0.0495
escan.add_potential         = "~/SiGe/pseudos/vq.Ge",      "~/SiGe/pseudos/vwr.pso", 0,0,0,0.2774
escan.add_potential         = "~/SiGe/pseudos/vq.SiGe.Si", "~/SiGe/pseudos/vwr.pso", 0,0,0,0.0495
escan.add_potential         = "~/SiGe/pseudos/vq.SiGe.Ge", "~/SiGe/pseudos/vwr.pso", 0,0,0,0.2774
escan.maskr                 = "~/SiGe/pseudos/maskr"
escan.references            = (-0.45, 0.10)
escan.vff = vff
escan.center    = (0,0,0) 

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
pools = 5
""" Number of pools of processes.

    During the evaluation of new individuals, the processes are splitted
    amongst this many pools. Each individual is computed by a single pool.
    If there are restrictions (as for pescan) as too many procs can be used by
    each calculation, they  are for the user to figure out.
"""

scales = [5.45, 5.55, 5.65]
""" Substrate lattice parameter. """
trials = range(3)
""" Number of trials per configuration space. """
ranges = [(8, 32)]
""" Size of the period in the configuration space. """
growth_direction = (0,0,1)
""" Epitaxial direction of the superlattice. """
dosym = 1
""" How to handle translation/inversion symmetries by individual. """


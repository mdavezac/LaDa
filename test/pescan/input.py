# sets up input stuff
vff = Vff()
vff.lattice.set_types = ("Si", "Ge"), ("Si","Ge")
vff.lattice.scale = 5.45
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

escan = Escan()
escan.eref                  = None
escan.cutoff                = 8.2
escan.smooth                = 1.0
escan.kinetic_scaling       = 1.0
escan.nbstates              = 10
escan.itermax               = 10
escan.nlines                = 50
escan.tolerance             = 1e-12
escan.rspace_cutoff         = 5e0
escan.fft_mesh              = 4*20, 20, int(ceil(sqrt(0.5) * 20))
escan.do_genpot             = True
escan.do_relax              = True
escan.input_wavefunctions   = None
escan.kpoint                = array([0,0,0])
escan.dnc_mesh              = None
escan.overlap_mesh          = None
escan.potential             = soH
escan.add_potential         = "test_input/pseudos/vq.Si", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.08
escan.add_potential         = "test_input/pseudos/vq.Ge", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.2774
escan.add_potential         = "test_input/pseudos/vq.SiGe.Si", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.08
escan.add_potential         = "test_input/pseudos/vq.SiGe.Ge", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.2774
escan.maskr                 = "test_input/pseudos/maskr"
escan.vff = vff

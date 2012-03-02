from quantities import angstrom

structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43)\
                     .add_atom(0,0,0, "Si")\
                     .add_atom(0.25,0.25,0.25, "Si")


vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
# vasp.lorbit     = 10
# vasp.npar       = 2
# vasp.lplane     = True
# vasp.addgrid    = True
# vasp.restart_from_contcar = False
vasp.set_smearing   = "metal", 0.01
vasp.relaxation = "volume", 50, 2
# vasp.set_symmetries = "on"

vasp.add_specie = "Si", "pseudos/Si"

# first_trial = { "encut": 0.9 }
# relaxation_dof = "volume ionic cellshape"
# relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)

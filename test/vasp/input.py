from quantities import angstrom
from lada.crystal.binary import zinc_blende

lattice = zinc_blende()
for site in lattice.sites: site.type = "Si"
lattice.scale = 5.43


vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n4 4 4\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1.5
vasp.set_smearing   = "metal", 0.01
vasp.relaxation = "volume", 50, 2
# vasp.set_symmetries = "on"

vasp.species = load("vasp/fere/pseudos", "Si", "Ge")

first_trial = { "encut": 0.9 }
relaxation_dof = "volume ionic cellshape"
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)

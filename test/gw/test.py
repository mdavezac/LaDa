from pylada.crystal.binary import zinc_blende
from pylada.mpi import world
from pylada.vasp import GWVasp

lattice = zinc_blende()
lattice.scale = 5.45
for site in lattice.sites: site.type = 'Si'

functional = GWVasp()
""" VASP functional """
functional.kpoints      = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
functional.precision    = "accurate"
functional.ediff        = 1e-5
functional.encut        = 0.8
functional.lorbit       = 10
functional.npar         = 2
functional.lplane       = True
functional.addgrid      = True
functional.set_smearing = "metal", 0.01
functional.relaxation   = "static"
functional.nbands       = 20
functional.vasp_library = "libvasp-5.2.11.so"

functional.add_specie = "Si", "pseudos/Si"

result = functional(lattice.to_structure(), outdir="results", comm=world)

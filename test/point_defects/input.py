""" NLEP fitting input script. """
from os.path import join, abspath
import spinel


workdir = None  #join(environ["SCRATCH"], environ["PBS_JOBNAME"])

lattice = spinel.lattice()
""" Sets back-bone lattice. """

supercell = array([[1, 0, 0],\
                   [0, 1, 0],\
                   [0, 0, 1]], dtype="float64" )
""" Supercell of defect structures. """

vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\ngamma\n6 6 10\n0 0 0"
vasp.precision  = "accurate"
vasp.smearing   = "bloechl"
vasp.ediff      = 1e-5
vasp.relaxation = "ionic"
vasp.encut      = 1.3
vasp.species    = species
vasp.workdir    = workdir
vasp.lorbit     = 10
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True

#                Symbol, directory of POTCAR, U parameters, max/min oxidation state, is magnetic
vasp.add_specie = "Rh", "pseudos/Rh", U=U("liechtenstein", "d", 3.3), 3, True
vasp.add_specie = "Zn", "pseudos/Zn", U=U("liechtenstein", "d", 6.0), 2, False
vasp.add_specie =  "O",  "pseudos/O", None, -2, False

# name of root output directories
outdir = "%s2%s%s4" % tuple([u for u in vasp.species.keys()])

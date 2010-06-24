""" NLEP fitting input script. """
from os import getcwd
from os.path import join, abspath, relpath
from lada.opt.changedir import Changedir
import spinel

outdir = abspath("results")
""" Root output directory. """
workdir = join(outdir, "tempdir")
""" Temporary calculation directory. """
if "NERSC_HOST" in environ:
  if environ["NERSC_HOST"] == "hopper":
    with Changedir(environ["HOME"]) as cwd:
      workdir = join(environ["SCRATCH"], relpath(outdir, getcwd()))

lattice = spinel.lattice()
""" Back-bone lattice. """
# changes species in lattice.
for site in lattice.sites:
  if   "A" in site.type: site.type[0] = "Rh"
  elif "B" in site.type: site.type[0] = "Zn"
  elif "X" in site.type: site.type[0] =  "O"

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
vasp.workdir    = workdir
vasp.lorbit     = 10
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
 
#                Symbol, directory of POTCAR, U parameters, max/min oxidation state, is magnetic
vasp.add_specie = "Rh", "pseudos/Rh", U("liechtenstein", "d", 3.3), 3, True
vasp.add_specie = "Zn", "pseudos/Zn", U("liechtenstein", "d", 6.0), 2, False
vasp.add_specie =  "O",  "pseudos/O", None, -2, False


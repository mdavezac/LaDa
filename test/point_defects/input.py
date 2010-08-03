""" NLEP fitting input script. """
from os import getcwd
from os.path import join, abspath, relpath
from lada.opt.changedir import Changedir
import spinel

lattice = spinel.lattice()
""" Back-bone lattice. """
# changes species in lattice.
for site in lattice.sites:
  if   "A" in site.type: site.type[0] = "Rh"
  elif "B" in site.type: site.type[0] = "Zn"
  elif "X" in site.type: site.type[0] =  "O"
lattice.scale = 8.506

supercell = array([[1, 0, 0],\
                   [0, 1, 0],\
                   [0, 0, 1]], dtype="float64" )
""" Supercell of defect structures. """

vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
vasp.lorbit     = 10
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.set_smearing   = "metal", 0.01
vasp.set_relaxation = "ionic"
vasp.set_symmetries = "off"

#                Symbol, directory of POTCAR, U parameters, max/min oxidation state, is magnetic
vasp.add_specie = "Rh", "pseudos/Rh", U("liechtenstein", "d", 3.3), 3, True
vasp.add_specie = "Zn", "pseudos/Zn", U("liechtenstein", "d", 6.0), 2, False
vasp.add_specie =  "O",  "pseudos/O", None, -2, False

queue     = "regular"
""" PBS queue. """
mppwidth  = 56
""" Total number of processes (eg core*nodes). """
walltime = "05:45:00"
""" PBS walltime. """
pools =  1
""" Number of process pools over which to parallelize each pbs scripts. 

    In other words, each vasp calculation will be launched on mppwidth/procpools cores.
"""
outdir    = "Zn2RhO4"
""" Root directory where to store results. """
relative  = "SCRATCH"
""" Relative calculation directory. 

    If not None, will use perform calculation in relative/directory, where
    directory is defined such that HOME/dictory is the outdir, and relative is
    an environment variable.
    In other words, if C{relative="SCRATCH"} and C{os.path.abspath(outdir) =
    $HOME/here/there}, then calculations will be performed in
    C{$SCRATCH/here/there}.
"""

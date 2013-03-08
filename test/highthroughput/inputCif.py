from glob import iglob
from pylada.crystal.read import icsd_cif
print "  test/hi/inputCif: entry"
vasp = Relax()
print "  test/hi/inputCif: === vasp ===\n%s\n=== end vasp ===" % (vasp,)


vasp.prec       = "accurate"
vasp.ediff      = 6.0e-5        # total, not per atom
vasp.encut      = 340.0
vasp.npar       = 1
vasp.ncore      = 1
vasp.lplane     = True

vasp.addgrid    = True
vasp.ismear     = "gaussian"
vasp.sigma      = 0.01
vasp.isym       = 0
vasp.relaxation = "ionic"
vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n20"
vasp.lorbit         = 10
vasp.add_param      = "lmaxmix",4

vasp.lcharg = True
vasp.lwave = True
vasp.lmaxmix = 4
vasp.loptics = False
vasp.lpead = False


# Fast: for fast, low accuracy use:
# vasp.prec           = "low"
# vasp.ediff          = 1e-3
# vasp.encut          = 0.5
# vasp.kpoints        = "\n0\nGamma\n1 1 1\n0. 0. 0.\n"
# Or could try:
# vasp.kpoints        = "\n0\nAuto\n20\n"


# See vasp/functional.py:  elementName, fileName, max or min oxidation state

pseudoDir = "/scratch/vnsteva/lada/pseudos"

vasp.add_specie = "Ag", pseudoDir + "/Ag"
vasp.add_specie = "Al", pseudoDir + "/Al"
vasp.add_specie = "As", pseudoDir + "/As"
vasp.add_specie = "Au", pseudoDir + "/Au"
vasp.add_specie = "B", pseudoDir + "/B"
vasp.add_specie = "Ba", pseudoDir + "/Ba"
vasp.add_specie = "Be", pseudoDir + "/Be"
vasp.add_specie = "Bi", pseudoDir + "/Bi"
vasp.add_specie = "Br", pseudoDir + "/Br"
vasp.add_specie = "C", pseudoDir + "/C"
vasp.add_specie = "Ca", pseudoDir + "/Ca"
vasp.add_specie = "Cd", pseudoDir + "/Cd"
vasp.add_specie = "Cl", pseudoDir + "/Cl"
vasp.add_specie = "Co", pseudoDir + "/Co"
vasp.add_specie = "Cr", pseudoDir + "/Cr"
vasp.add_specie = "Cu", pseudoDir + "/Cu"
vasp.add_specie = "F", pseudoDir + "/F"
vasp.add_specie = "Fe", pseudoDir + "/Fe"
vasp.add_specie = "Ga", pseudoDir + "/Ga"
vasp.add_specie = "Ge", pseudoDir + "/Ge"
vasp.add_specie = "H", pseudoDir + "/H"
vasp.add_specie = "Hf", pseudoDir + "/Hf"
vasp.add_specie = "Hg", pseudoDir + "/Hg"
vasp.add_specie = "I", pseudoDir + "/I"
vasp.add_specie = "In", pseudoDir + "/In"
vasp.add_specie = "Ir", pseudoDir + "/Ir"
vasp.add_specie = "K", pseudoDir + "/K"
vasp.add_specie = "La", pseudoDir + "/La"
vasp.add_specie = "Li", pseudoDir + "/Li"
vasp.add_specie = "Lu", pseudoDir + "/Lu"
vasp.add_specie = "Mg", pseudoDir + "/Mg"
vasp.add_specie = "Mn", pseudoDir + "/Mn"
vasp.add_specie = "Mo", pseudoDir + "/Mo", U("dudarev", "d", 3.0)
vasp.add_specie = "N", pseudoDir + "/N"
vasp.add_specie = "Na", pseudoDir + "/Na"
vasp.add_specie = "Nb", pseudoDir + "/Nb"
vasp.add_specie = "Ni", pseudoDir + "/Ni"
vasp.add_specie = "O", pseudoDir + "/O"
vasp.add_specie = "P", pseudoDir + "/P"
vasp.add_specie = "Pb", pseudoDir + "/Pb"
vasp.add_specie = "Pd", pseudoDir + "/Pd"
vasp.add_specie = "Pt", pseudoDir + "/Pt"
vasp.add_specie = "Rb", pseudoDir + "/Rb"
vasp.add_specie = "Rh", pseudoDir + "/Rh"
vasp.add_specie = "Ru", pseudoDir + "/Ru"
vasp.add_specie = "S", pseudoDir + "/S"
vasp.add_specie = "Sb", pseudoDir + "/Sb"
vasp.add_specie = "Sc", pseudoDir + "/Sc"
vasp.add_specie = "Se", pseudoDir + "/Se"
vasp.add_specie = "Si", pseudoDir + "/Si"
vasp.add_specie = "Sn", pseudoDir + "/Sn"
vasp.add_specie = "Sr", pseudoDir + "/Sr"
vasp.add_specie = "Ta", pseudoDir + "/Ta"
vasp.add_specie = "Te", pseudoDir + "/Te"
vasp.add_specie = "Ti", pseudoDir + "/Ti"
vasp.add_specie = "V", pseudoDir + "/V"
vasp.add_specie = "W", pseudoDir + "/W"
vasp.add_specie = "Y", pseudoDir + "/Y"
vasp.add_specie = "Zn", pseudoDir + "/Zn"
vasp.add_specie = "Zr", pseudoDir + "/Zr"


vasp.species["Mo"].moment = [1.e0,4.e0]

vasp.first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }
""" parameter to override during first relaxation step. """

vasp.relaxation = "volume ionic cellshape"
""" Degrees of freedom to relax. """

vasp.maxiter = 5
""" Maximum number of iterations before bailing out. """

vasp.keep_steps = True
""" Whether to keep or delete intermediate steps. """


matLatPairs = []

for fname in iglob("icsd_structures/*.cif"):
  material = fname[fname.index('/')+1:-4]
  lattice = icsd_cif( fname)
  print "  test/hi/inputCif: material: ", material
  print "    test/hi/inputCif: lattice: ", lattice
  matLatPairs.append( (material, lattice))

print "test/hi/inputFixedCif: matLatPairs len: %d" % (len( matLatPairs),)




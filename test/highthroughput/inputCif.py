from glob import iglob
from pylada.crystal.read import icsd_cif
print "  test/hi/inputCif: entry"
vasp = Relax()
print "  test/hi/inputCif: === vasp ===\n%s\n=== end vasp ===" % (vasp,)

vasp.precision      = "accurate"
vasp.ediff          = 1e-5 # precision per ATOM
vasp.encut          = 1.0


# Fast: change precision from accurate to low
vasp.precision      = "low"

# Fast: change ediff from 1.e-5 to 1.e-3
vasp.ediff          = 1e-3 # precision per ATOM

# Fast: change encut from 1.0 to 0.5
vasp.encut          = 0.5

# Fast: specify kpoints
vasp.kpoints        = "\n0\nGamma\n1 1 1\n0. 0. 0.\n"
# Or could try:
#vasp.kpoints        = "\n0\nAuto\n10\n"


# See vasp/functional.py:  elementName, fileName, max or min oxidation state
vasp.add_specie = "Al", "pseudos/Al", None, 3
vasp.add_specie = "Mg", "pseudos/Mg", None, 2
vasp.add_specie = "Mo", "pseudos/Mo"
vasp.add_specie = "O", "pseudos/O", None, -2
vasp.add_specie = "S", "pseudos/S", None, -2

#first_trial = { "kpoints": "\n0\nAuto\n40", "encut": 1.0 }
vasp.first_trial = {}
""" parameter to override during first relaxation step. """
vasp.relaxation = "volume ionic cellshape"
""" Degrees of freedom to relax. """
vasp.maxiter = 5
""" Maximum number of iterations before bailing out. """
vasp.keep_steps = True
""" Whether to keep or delete intermediate steps. """


def scale(structure):
  """ Returns *guessed* scale (eg volume^(0.33)) for a given structure. """
  from numpy.linalg import det
  if "O" in [atom.type for atom in structure]:    spvol = 8.5**3/4e0
  elif "Se" in [atom.type for atom in structure]: spvol = 9.5**3/4e0
  elif "Te" in [atom.type for atom in structure]: spvol = 10.5**3/4e0
  else: raise ValueError("unknown atom.type: %s" % (atom.type,))

  nfu = float(len(structure)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0)


matLatPairs = []

for fname in iglob("icsd_structures/*.cif"):
  material = fname[fname.index('/')+1:-4]
  lattice = icsd_cif( fname)
  print "  test/hi/inputCif: material: ", material
  print "    test/hi/inputCif: lattice: ", lattice
  matLatPairs.append( (material, lattice))

print "test/hi/inputFixedCif: matLatPairs len: %d" % (len( matLatPairs),)




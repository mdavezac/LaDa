from pylada.crystal import A2BX4
print "  test/hi/input: entry"
vasp = Relax()
print "  test/hi/input: === vasp ===\n%s\n=== end vasp ===" % (vasp,)

vasp.precision      = "accurate"
vasp.ediff          = 1e-5 # precision per ATOM
vasp.encut          = 1.0

vasp.add_specie = "Al", "pseudos/Al", None, 3
vasp.add_specie = "Mg", "pseudos/Mg", None, 2
vasp.add_specie = "O", "pseudos/O", None, -2

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

  nfu = float(len(structure)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0)

#problem with Rh2TiO4 for some reason

""" Materials to compute. """
materials = [ "Al2MgO4"]

""" Number of random anti-ferro trials. """
nbantiferro = 8
nbrandom    = 3
do_ferro    = False
do_antiferro = False

lattices = [A2BX4.b5(), A2BX4.b21()]

mlen = len( materials)
llen = len( lattices)
matLatPairs = (mlen * llen) * [None]
print "  test/hi/input: mlen: ", mlen
print "  test/hi/input: llen: ", llen
print "  test/hi/input: pairs len: ", len(matLatPairs)
kk = 0
for mat in materials:
  print "  test/hi/input: mat: ", mat
  for lat in lattices:
    print "    test/hi/input: lat: ", lat
    matLatPairs[kk] = (mat, lat,)
    kk += 1

print "test/hi/inputFixed: mats len: %d  lats len: %d  matLatPairs len: %d" \
  % (mlen, llen, len( matLatPairs),)


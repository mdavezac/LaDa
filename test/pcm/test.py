from lada.pcm import Clj, bond_name
from numpy import array

from lada.opt import read_input
from lada.crystal import A2BX4
from lada.crystal import fill_structure, Structure

clj  = Clj()
""" Point charge + r^12 + r^6 model. """
clj.ewald_cutoff = 32

clj.charges["A"] = 3.0
clj.charges["B"] = 2.0
clj.charges["X"] = -2.0

# clj.lj_cutoff = 10.5
# clj.mesh = (15, 15, 15)
# hs = { "A":1.34, "B":2.25,"X":1.14}
# vdw = {"A":2.2,"B":3,"X":1.9}
# for a in ["A", "B", "X" ]:
#   for b in ["A", "B", "X" ]:
#     type = bond_name( a, b )
#     hs_ = float( hs[a] ) + float( hs[b]  )
#     vdw_ = float(vdw[a]) + float(vdw[b] )
#     clj.bonds[type] = pow(hs_, 12.0), pow(vdw_, 6.0)



x = 0.387
lattice = A2BX4.b5(x)
lattice.scale = 8.5
structure = fill_structure(lattice.cell, lattice)

N = clj.ewald(structure)
nbB = 2
for atom in structure.atoms:
  if atom.type == "X": continue
  elif atom.type == "A" and nbB > 0:
    atom.type = "B"
    nbB -= 1
  elif atom.type == "B": atom.type = "A"

I = clj.ewald(structure)
print I.energy, N.energy
print I.energy - N.energy,  structure.scale, x

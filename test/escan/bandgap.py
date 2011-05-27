from os.path import join
from numpy import matrix
from lada.escan import BandGap, read_input
from lada.crystal import FreezeCell

# reads input file.
input = read_input("input.py")

# creating unrelaxed structure.
cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
structure = input.vff.lattice.to_structure(cell)
for i, atom in enumerate(structure.atoms):
  atom.type = "Si" if i % 10 < 6 else "Ge"
structure.scale = 5.65
structure.freeze = FreezeCell.a0 | FreezeCell.a1

functional = BandGap(escan=input.escan, references=(0.2, -0.5))
out = functional(structure, outdir=join("results", "BG"))
print "%f - %f = %f " % (out.cbm, out.vbm, out.bandgap)

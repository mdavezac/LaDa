from lada.crystal import Structure
from lada.pcm import Clj, bond_name
from lada.physics import a0, Ry
from quantities import angstrom, eV, hartree

clj  = Clj()
""" Point charge + r^12 + r^6 model. """
clj.ewald_cutoff = 100

clj.charges["A"] = -1.0
clj.charges["B"] =  1.0

structure = Structure()
structure.set_cell = (1,0,0),\
                     (0,1,0),\
                     (0,0,1)
structure.scale = 2
structure.add_atom = (0,0,0), "A"
structure.add_atom = (a0.rescale(angstrom)/structure.scale,0,0), "B"

print clj.ewald(structure).energy, hartree.rescale(eV)

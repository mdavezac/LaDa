from lada.crystal import binary, supercell
from lada.vff import Functional
a = binary.zinc_blende()
a = supercell(binary.zinc_blende(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]])
b = Functional.build_tree(a, overlap=0.5)

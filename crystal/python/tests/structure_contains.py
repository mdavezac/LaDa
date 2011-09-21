from numpy import all, abs, array, any
from collections import namedtuple
from lada.crystal import Structure
from lada.crystal._crystal import StructureStr, StructureVec, StructureSet

Atom = namedtuple('Atom', ['type', 'pos'])


cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, weight=0e0)
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Si")
assert "Si" in a
assert "Ge" not in a
assert (0.25, 0.25, 0.25) in a
assert (0.75, 0.75, 0.25) in a
assert (0.7, 0.75, 0.25) not in a
assert Atom(type="Si", pos=(-0.75, 0.25, -0.75)) in a
assert Atom(type="Ge", pos=(-0.75, 0.25, -0.75)) not in a
assert Atom(type="Si", pos=(-0.7, 0.25, -0.75)) not in a

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, weight=0e0, kind="list")
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Si", "Ge")
assert ["Si"] in a
assert ["Ge"] not in a
assert ["Ge", "Si"] in a
assert [] not in a
assert (0.25, 0.25, 0.25) in a
assert (0.75, 0.75, 0.25) in a
assert Atom(type=["Si"], pos=(-0.75, 0.25, -0.75)) not in a
assert Atom(type=["Si", "Ge"], pos=(-0.75, 0.25, -0.75)) in a
assert Atom(type=["Si", "Ge"], pos=(-0.7, 0.25, -0.75)) not in a
a.add_atom(0,0,0,[])
assert [] in a

cell = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
a = Structure(cell, weight=0e0, kind="set")
a.add_atom(0,0,0, "Si")\
          (0.25,0.25,0.25, "Si", "Ge")
assert ["Si"] in a
assert ["Ge"] not in a
assert ["Ge", "Si"] in a
assert (0.25, 0.25, 0.25) in a
assert (0.75, 0.75, 0.25) in a
assert Atom(type=[], pos=(-0.75, 0.25, -0.75)) not in a
assert Atom(type=["Si"], pos=(-0.75, 0.25, -0.75)) not in a
assert Atom(type=["Si", "Ge"], pos=(-0.75, 0.25, -0.75)) in a
assert Atom(type=["Si", "Ge"], pos=(-0.7, 0.25, -0.75)) not in a
a.add_atom(0,0,0,[])
assert Atom(type=[], pos=(0, 0, 0)) in a

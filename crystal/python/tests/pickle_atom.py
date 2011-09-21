from pickle import loads, dumps
from numpy.linalg import norm
from lada.crystal import Atom
from lada.crystal._crystal import AtomStr, AtomVec, AtomSet

a = Atom(0.25,0.25,0,"Au", site=-1, freeze=1)
string = dumps(a)
b = loads(string)
assert type(b) == AtomStr
assert norm(a.pos-b.pos) < 1e-12
assert a.site == b.site
assert a.type == b.type
assert a.freeze == b.freeze

a = Atom(0.25,0.25,0,["Au", "Pd"], site=-1, freeze=1)
string = dumps(a)
b = loads(string)
assert type(b) == AtomVec
assert norm(a.pos-b.pos) < 1e-12
assert a.site == b.site
assert a.type == b.type
assert a.freeze == b.freeze

a = Atom(0.25,0.25,0,set(["Au", "Pd"]), site=-1, freeze=1)
string = dumps(a)
b = loads(string)
assert type(b) == AtomSet
assert norm(a.pos-b.pos) < 1e-12
assert a.site == b.site
assert a.type == b.type
assert a.freeze == b.freeze

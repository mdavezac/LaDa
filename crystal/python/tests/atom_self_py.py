# Check that wrappers around a c++-created instance work and that the dictionary is
# conserved for the life of the c++ object, rather that of the python wrapper.
from numpy import abs, all
from lada.crystal.cppwrappers.atom import @PYTYPE@
from atom_self_@TYPE@ import get_new_object, get_static_object
a = get_static_object()
a.pos = [0.1, 0.1, 0.1] # set the whole array.
a.pos[1] = 0.2          # set a single item.
a.type = 'Au' if "@PYTYPE@" == "AtomStr" else ['Au', 'Pd'] # set whole array
a.m = 0.5 # add other attribute.
i = id(a.__dict__) # keep track of __dict__'s id to check identity later on.
del a # delete wrapper.
b = get_static_object() # gets the same c++ object back
assert all(abs(b.pos - [0.1, 0.2, 0.1]) < 1e-12)
assert b.type == 'Au' if "@PYTYPE@" == "AtomStr" else ['Au', 'Ga']
assert id(b.__dict__) == i
assert hasattr(b, 'm')
assert abs(b.m - 0.5) < 1e-12


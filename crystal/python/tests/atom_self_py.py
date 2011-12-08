# Check that wrappers around a c++-created instance work and that the dictionary is
# conserved for the life of the c++ object, rather that of the python wrapper.
from numpy import abs, all
from lada.crystal.cppwrappers.atom import @PYTYPE@
from atom_self_@TYPE@ import get_new_object, get_static_object, self_is_null

# Returns wrapper to a static object in atom_self extension.
a = get_static_object()
# check static objects reference to wrapper is not null.
assert not self_is_null()
# Gets second wrapper.
b = get_static_object()
# check that the two wrappers are the same object. Only one living wrapper at a time.
assert a is b
del b
# checks deleting b does not invalidate reference to a. atomsomethin_dealloc
# should not have been called.
assert not self_is_null()

# changes the atom, including a python specific attribute.
a.pos = [0.1, 0.1, 0.1] # set the whole array.
a.pos[1] = 0.2          # set a single item.
a.type = 'Au' if "@PYTYPE@" == "AtomStr" else ['Au', 'Pd'] # set whole array
a.m = 0.5 # add other attribute.
i = id(a.__dict__) # keep track of __dict__'s id to check identity later on.
del a # delete wrapper.
# There should be no reference left to wrapper inside the static object.
assert self_is_null()
# get another wrapper around static object. This time it is new.
b = get_static_object()
# check reference inside static object has been set again.
assert not self_is_null()
# makes all values, including the python specific attribute, have been retained.
assert all(abs(b.pos - [0.1, 0.2, 0.1]) < 1e-12)
assert b.type == 'Au' if "@PYTYPE@" == "AtomStr" else ['Au', 'Ga']
assert id(b.__dict__) == i
assert hasattr(b, 'm')
assert abs(b.m - 0.5) < 1e-12
del b
assert self_is_null()


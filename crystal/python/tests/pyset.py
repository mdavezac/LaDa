from lada.crystal.cppwrappers.atom import _new_set as SpecieSet

# checks initialization and some equalities.
a = SpecieSet(['Au', 'Pd'])
# assert a == set()
assert set(['Au']) != a
print "HERE"
# assert a != set(['Au', 'Pd'])
# assert set(['Au', 'Pd']) != a

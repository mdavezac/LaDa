from lada.crystal.cppwrappers.atom import _new_set as SpecieSet

# checks initialization and some equalities.
a = SpecieSet(['Au', 'Pd'])
class NewSet(set):
  def __eq__(self, a): 
    print "IN DERIVED", self, super(NewSet, self).__eq__(set(['Au', 'Pd']))
    return super(NewSet, self).__eq__(a)
class NewSetB(set):
  def __eq__(self, a): 
    print "IN DERIVED B", self, super(NewSet, self).__eq__(set(['Au', 'Pd']))
    return super(NewSet, self).__eq__(a)
# assert a == set()
assert set(['Au', 'Pd']) == a
print "HERE"
assert NewSet(['Au', 'Pd']) == NewSetB(['Au', 'Pd'])
print "THERE"
# assert a != set(['Au', 'Pd'])
# assert set(['Au', 'Pd']) != a

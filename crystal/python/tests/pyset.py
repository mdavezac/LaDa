from lada.crystal.cppwrappers.atom import _new_set as SpecieSet

# checks initialization and some equalities.
a = SpecieSet()
assert a == set()
assert a != set(['Au', 'Pd'])
a = SpecieSet(["Au", "Pd"][i] for i in range(2))
assert a == set(['Au', 'Pd'])
assert a != set(['Au', 'Pd', 'Ga'])
a = SpecieSet(["Au", "Pd"])
assert a == SpecieSet(['Au', 'Pd'])
assert a != SpecieSet(['Au', 'Pd', 'Ga'])

assert SpecieSet(["Au", "Pd"]) != []
assert SpecieSet(["Au", "Pd"]) == ['Au', 'Pd']
assert SpecieSet([]) == []
assert [] == SpecieSet([])
assert [] != SpecieSet(["Au", "Pd"])
assert ['Au', 'Pd'] == SpecieSet(["Au", "Pd"])
assert SpecieSet(["Au", "Pd"]) != ["Au"]
assert SpecieSet(["Au", "Pd"]).isdisjoint(["Au"]) == False
assert SpecieSet(["Au", "Pd"]).isdisjoint(["Cr"]) == True
assert SpecieSet(["Au", "Pd"]).issubset(["Cr", "Ag"]) == False
assert SpecieSet(["Au", "Pd"]).issubset(["Au"]) == False
assert SpecieSet(["Au", "Pd"]).issubset(["Au", "Pd"]) == True
assert SpecieSet(["Au", "Pd"]).issubset(["Au", "Au", "Pd"]) == True
assert SpecieSet(["Au", "Pd"]).issubset(["Au", "Pd", "Cr"]) == True
assert (SpecieSet(["Au", "Pd"]) <= ["Cr", "Ag"]) == False
assert (SpecieSet(["Au", "Pd"]) <= ["Au"]) == False
assert (SpecieSet(["Au", "Pd"]) <= ["Au", "Pd"]) == True
assert (SpecieSet(["Au", "Pd"]) <= ["Au", "Au", "Pd"]) == True
assert (SpecieSet(["Au", "Pd"]) <= ["Au", "Pd", "Cr"]) == True
assert (SpecieSet(["Au", "Pd"]) < ["Au", "Pd", "Cr"]) == True
assert (SpecieSet(["Au", "Pd"]) < ["Au", "Pd"]) == False
assert (SpecieSet(["Au", "Pd"]) < ["Au", "Au"]) == False
assert (SpecieSet(["Au", "Pd"]) > ["Au", "Pd", "Cr"]) == False
assert (SpecieSet(["Au", "Pd"]) > ["Au", "Pd"]) == False
assert (SpecieSet(["Au", "Pd"]) > ["Au", "Au"]) == False
assert (SpecieSet(["Au", "Pd"]) > ["Au"]) == True
assert (SpecieSet(["Au", "Pd"]) >= ["Au"]) == True
assert (SpecieSet(["Au", "Pd"]) >= ["Cr"]) == False
assert (SpecieSet(["Au", "Pd"]) >= ["Au", "Pd", "Cr"]) == False

assert (SpecieSet(["Au", "Pd"]) | ["Cr"]) == set(["Au", "Pd", "Cr"])
assert (SpecieSet(["Au", "Pd"]) | ["Cr", "Ag"]) != set(["Au", "Pd", "Cr"])
assert (SpecieSet(["Au", "Pd"]) | []) == set(["Au", "Pd"])
a = SpecieSet(["Au", "Pd"])
b = a
a |= ["Cr"] 
assert b == set(["Au", "Pd", "Cr"])
a = SpecieSet(["Au", "Pd"])
b = a
a |= []
assert b == set(["Au", "Pd"])

assert SpecieSet() == []
assert (SpecieSet(["Au", "Pd"]) & ["Cr"]) == set([])
assert (SpecieSet(["Au", "Pd"]) & ["Pd"]) == set(["Pd"])
a = SpecieSet(["Au", "Pd"])
b = a
a &= set(["Cr"])
assert b == []
a = SpecieSet(["Au", "Pd"])
b = a
a &= ["Pd"] 
assert b == ["Pd"]

assert SpecieSet(["Au", "Pd"]) - ["Au", "Cr"] == set(["Pd"])
a = SpecieSet(["Au", "Pd"])
b = a
a -= ["Au", "Cr"] 
assert b == ["Pd"]

assert SpecieSet(["Au", "Pd"]) ^ ["Au", "Cr"] == set(["Pd", "Cr"])
a = SpecieSet(["Au", "Pd"])
b = a
a ^= ["Au", "Cr"] 
assert b == ["Pd", "Cr"]

a = SpecieSet(["Au", "Pd"])
a.add("Au")
assert a == ["Au", "Pd"]
a.add("Cr")
assert a == ["Au", "Pd", "Cr"]
a.remove("Cr")
assert a == ["Au", "Pd"]
try: a.remove("Cr")
except KeyError: pass
try: a.remove("Cr")
except IndexError: pass

assert a.pop("Au") == "Au"
a.clear()
assert a == []
assert str(a) == "set()"
assert str(SpecieSet(["Au", "Pd"])) == "{\"Au\", \"Pd\"}"

a = SpecieSet(["Au", "Pd"])
a.discard("Cr")
a.discard("Au")
assert a == ["Pd"]

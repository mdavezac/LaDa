from lada.crystal._crystal import SpecieSet
from lada.error import out_of_range

a = SpecieSet()
a = SpecieSet(["Au", "Pd"][i] for i in range(2))
a = SpecieSet(["Au", "Pd"])

assert SpecieSet(["Au", "Pd"]) != []
assert SpecieSet([]) == []
assert [] == SpecieSet([])
assert [] != SpecieSet(["Au", "Pd"])
assert SpecieSet(["Au", "Pd"]) != ["Au"]
assert SpecieSet(["Au", "Pd"]) != ["Au", "Pd"]
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

assert (SpecieSet(["Au", "Pd"]) | ["Cr"]) == ["Au", "Pd", "Cr"]
assert (SpecieSet(["Au", "Pd"]) | ["Cr", "Ag"]) != ["Au", "Pd", "Cr"]
assert (SpecieSet(["Au", "Pd"]) | []) == ["Au", "Pd"]
a = SpecieSet(["Au", "Pd"])
b = a
a |= ["Cr"] 
assert b == ["Au", "Pd", "Cr"]
a = SpecieSet(["Au", "Pd"])
b = a
a |= []
assert b == ["Au", "Pd"]

assert SpecieSet() == []
assert (SpecieSet(["Au", "Pd"]) & ["Cr"]) == []
assert (SpecieSet(["Au", "Pd"]) & ["Pd"]) == ["Pd"]
a = SpecieSet(["Au", "Pd"])
b = a
a &= ["Cr"] 
assert b == []
a = SpecieSet(["Au", "Pd"])
b = a
a &= ["Pd"] 
assert b == ["Pd"]

assert SpecieSet(["Au", "Pd"]) - ["Au", "Cr"] == ["Pd"]
a = SpecieSet(["Au", "Pd"])
b = a
a -= ["Au", "Cr"] 
assert b == ["Pd"]

assert SpecieSet(["Au", "Pd"]) ^ ["Au", "Cr"] == ["Pd", "Cr"]
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
except out_of_range: pass
try: a.remove("Cr")
except IndexError: pass

assert a.pop("Au") == "Au"
a.clear()
assert a == []
assert str(a) == "SpecieSet()"
assert str(SpecieSet(["Au", "Pd"])) == "SpecieSet([\"Au\", \"Pd\"])"

a = SpecieSet(["Au", "Pd"])
a.discard("Cr")
a.discard("Au")
assert a == ["Pd"]

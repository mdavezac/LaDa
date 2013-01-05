from pylada.crystal.A2BX4 import b5
from pylada.enumeration import Enum
from spinel_check import checkme


# create spinel lattice and make it inverse.
lattice = Enum(b5())
for site in lattice.sites:
  if "X" in site.type: continue
  site.type = 'A' if 'B' in site.type else ('A', 'B')
lattice.find_space_group()

# look at inverse only.
def check_concentration(x, flavorbase, smith):
  from pylada.enumeration import as_numpy

  types = as_numpy(x, flavorbase) 
  result = 2*len(types[types == 1]) == len(types)

  return result, None if result else False

result = [i for i, dummy, dummy, dummy in lattice.xn(2)]
result.append([i for i, dummy, dummy, dummy in lattice.xn(4)])

# creates file with result.
# with open("spinel_check.py", "w") as file: 
#   file.write('checkme = {0}\n'.format(repr(result)))

for a, b in zip(checkme, result): assert a == b

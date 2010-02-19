from pickle import dumps, loads
from kpoints import Gamma, Density
from specie  import Specie
from launch  import Launch
def test_gamma_serialization():

  g = Gamma()

  string = dumps(g)
  e = loads(string)
  e.newattr = 0
  assert not hasattr(g, "newattr"), "Copied, not serialized."


def test_density_serialization():

  g = Density(offset=(0.4,0.5,0.5), length=20)
  string = dumps(g)
  g.offset = (0,0,0)
  e = loads(string)
  assert e.offset == (0.4,0.5,0.5) and e.length == 20

def test_species_serialization():
  from os.path import expanduser

  g = Specie("Ag", "~/usr/pseudos/Ag")

  string = dumps(g)
  g.symbol = "Au"
  e = loads(string)

  assert e.symbol == "Ag"
  assert e.path == expanduser("~/usr/pseudos/Ag")

def test_launch_serialization():
  from os.path import expanduser

  g = Specie("Ag", "~/usr/pseudos/Ag")

  string = dumps(g)
  g.symbol = "Au"
  e = loads(string)

  assert e.symbol == "Ag"
  assert e.path == expanduser("~/usr/pseudos/Ag")

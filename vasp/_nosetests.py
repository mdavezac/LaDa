from pickle import dumps, loads
from kpoints import Gamma, Density
from specie  import Specie, nlep, U
from launch  import Launch
from incar  import Standard, Incar
from .  import Vasp
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

  species = Specie\
            (
              "Sn", 
              path="~/nlep/pseudos/Sn", 
              U=[nlep(type="Dudarev", l=i, U0=0e0) for i in ["s", "p", "d"]]
            ),\
            Specie\
            (
              "O", 
              path="~/nlep/pseudos/O",
              U=[nlep(type="Dudarev", l=i, U0=0e0) for i in ["s", "p"]]
            )
  species[0].U.append( U(type="Dudarev", U=4e0, l=2) ) # add U to Sn atoms.

  string = dumps(species)
  e = loads(string)
  species[0].symbol = "Au"

  assert e[0].symbol == "Sn"
  assert e[0].path == expanduser("~/nlep/pseudos/Sn")

def test_incar_serialization():
  incar = Incar()
  string = dumps(incar)
  e = loads(string)
  incar.iniwave.value = "jellium"
  assert e.iniwave.value == 1, str(e.iniwave.value)

def test_launch_serialization():
  incar = Launch()
  string = dumps(incar)
  e = loads(string)
  incar.iniwave.value = "jellium"
  assert e.iniwave.value == 1, str(e.iniwave.value)

def test_vasp_serialization():
  from specie  import Specie, nlep, U
  species = Specie\
            (
              "Sn", 
              path="~/nlep/pseudos/Sn", 
              U=[nlep(type="Dudarev", l=i, U0=0e0) for i in ["s", "p", "d"]]
            ),\
            Specie\
            (
              "O", 
              path="~/nlep/pseudos/O",
              U=[nlep(type="Dudarev", l=i, U0=0e0) for i in ["s", "p"]]
            )
  """ Parameters of atomic species. """
  species[0].U.append( U(type="Dudarev", U=4e0, l=2) ) # add U to Sn atoms.
  vasp = Vasp\
         (
           kpoints    = "Automatic generation\n0\ngamma\n6 6 10\n0 0 0",
           precision  = "accurate",
           smearing   = "bloechl",
           ediff      = 1e-5,
           relaxation = "ionic",
           encut      = 1.3, # uses ENMAX * 1, which is VASP default
           species    = species,
           workdir    = None
         )
  """ VASP functional """
  # adds some extra parameters to VASP functional.
  vasp.nbands     = Standard("NBANDS", 64)
  vasp.lorbit     = Standard("LORBIT", 10)
  vasp.npar       = Standard("NPAR", 2)
  vasp.lplane     = Standard("LPLANE", ".TRUE.")
  vasp.addgrid    = Standard("ADDGRID", ".TRUE.")

  string = dumps(vasp)
  e = loads(string)
  vasp.iniwave.value = "jellium"
  assert e.iniwave.value == 1, str(e.iniwave.value)

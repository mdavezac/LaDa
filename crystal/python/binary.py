""" Common binary lattices. """

def zinc_blende():
  """ Zinc-blende lattice. """
  from . import Lattice

  lattice = Lattice()
  lattice.name = "Zinc-Blende"
  lattice.scale = 1e0
  lattice.set_cell = (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0)
  lattice.add_sites = [(   0,   0,   0), 'A'], \
                      [(0.25,0.25,0.25), 'B']
  lattice.find_space_group()
  return lattice

def rock_salt():
  """ Rock-salt lattice. """
  from . import Lattice

  lattice = Lattice()
  lattice.name = "Rock-Salt"
  lattice.scale = 1e0
  lattice.set_cell = (1e0, 0, 0), (0, 1e0, 0), (0, 0, 1e0)
  lattice.add_sites = [(  0,  0,  0), 'A'], \
                      [(0.5,0.5,0.5), 'B']
  lattice.find_space_group()
  return lattice

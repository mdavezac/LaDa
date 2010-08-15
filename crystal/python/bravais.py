""" Defines Bravais Lattices. """

def fcc(type="A", scale=1e0):
  """ Creates an FCC lattice with a single site. """
  from . import Lattice

  lattice = Lattice()
  lattice.scale = scale
  lattice.set_cell = (0e0, 0.5, 0.5),\
                     (0.5, 0.0, 0.5),\
                     (0.5, 0.5, 0.0)
  lattice.name = "fcc"
  lattice.add_site = (0e0, 0e0, 0e0), type
  lattice.find_space_group()
  return lattice

def bcc(type="A", scale=1e0):
  """ Creates an FCC lattice with a single site. """
  from . import Lattice

  lattice = Lattice()
  lattice.scale = scale
  lattice.set_cell = (-0.5,  0.5,  0.5),\
                     ( 0.5, -0.5,  0.5),\
                     ( 0.5,  0.5, -0.5)
  lattice.name = "bcc"
  lattice.add_site = (0e0, 0e0, 0e0), type
  lattice.find_space_group()
  return lattice

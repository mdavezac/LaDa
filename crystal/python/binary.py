""" Holds different binary lattices. """

def zinc_blende():
  """ Defines a default Zinc-Blende lattice. """
  from ..crystal import Lattice, Site
  from numpy import array

  lattice = Lattice()
  lattice.set_cell = (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0)
  lattice.scale = 1e0
  lattice.add_site = (0,0,0), "A"
  lattice.add_site = (0.25,0.25,0.25), "B"
  lattice.name = "Zinc-Blende"
  return lattice

""" Defines Bravais Lattices. """
__docformat__ = "restructuredtext en"

def fcc(type="A", scale=1e0):
  """ Creates an FCC lattice with a single site. 
  
      :Parameters:
        type : str
          A string defining the atomic-specied of this crystal structure.
        scale : float, quantity.Unit
          Scale in units of Bhor radius. If given a float with a dimension,
          will convert it to atomic units.
      :return: a `lada.crystal.Lattice` object for this crystal type.
  """
  from . import Lattice
  from ..physics import a0

  lattice = Lattice()
  lattice.scale = float(scale.rescale(a0) if hasattr(scale, "rescale") else scale)
  lattice.set_cell = (0e0, 0.5, 0.5),\
                     (0.5, 0.0, 0.5),\
                     (0.5, 0.5, 0.0)
  lattice.name = "fcc"
  lattice.add_site = (0e0, 0e0, 0e0), type
  lattice.find_space_group()
  return lattice

def bcc(type="A", scale=1e0):
  """ Creates an FCC lattice with a single site. 
  
      :Parameters:
        type : str
          A string defining the atomic-specied of this crystal structure.
        scale : float, quantity.Unit
          Scale in units of Bhor radius. If given a float with a dimension,
          will convert it to atomic units.
      :return: a `lada.crystal.Lattice` object for this crystal type.
  """
  from . import Lattice
  from ..physics import a0

  lattice = Lattice()
  lattice.scale = float(scale.rescale(a0) if hasattr(scale, "rescale") else scale)
  lattice.set_cell = (-0.5,  0.5,  0.5),\
                     ( 0.5, -0.5,  0.5),\
                     ( 0.5,  0.5, -0.5)
  lattice.name = "bcc"
  lattice.add_site = (0e0, 0e0, 0e0), type
  lattice.find_space_group()
  return lattice

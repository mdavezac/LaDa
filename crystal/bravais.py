""" Defines Bravais Lattices. """
__docformat__ = "restructuredtext en"
def bcc():
  """ Creates a BCC lattice with a single site. """
  from lada.crystal import Structure
  return Structure( -0.5, 0.5, 0.5,\
                    0.5, -0.5, 0.5,\
                    0.5, 0.5, -0.5,\
                    scale=1, name='bcc' )\
           .add_atom(0, 0, 0, 'A')

def fcc():
  """ Creates an FCC lattice with a single site. """
  from lada.crystal import Structure
  return Structure( 0, 0.5, 0.5,\
                    0.5, 0, 0.5,\
                    0.5, 0.5, 0,\
                    scale=1, name='fcc' )\
           .add_atom(0, 0, 0, 'A')

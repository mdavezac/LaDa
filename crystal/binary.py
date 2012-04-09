""" Defines binary lattices. """
__docformat__ = "restructuredtext en"

def rock_salt():
  """ rock_salt lattice """
  from lada.crystal import Structure
  return Structure( 1, 0, 0,\
                    0, 1, 0,\
                    0, 0, 1,\
                    scale=1, name='Rock-Salt' )\
           .add_atom(0, 0, 0, 'A')\
           .add_atom(0.5, 0.5, 0.5, 'B')

def zinc_blende():
  """ zinc_blende lattice """
  from lada.crystal import Structure
  return Structure( 0, 0.5, 0.5,\
                    0.5, 0, 0.5,\
                    0.5, 0.5, 0,\
                    scale=1, name='Zinc-Blende' )\
           .add_atom(0, 0, 0, 'A')\
           .add_atom(0.25, 0.25, 0.25, 'B')

def wurtzite():
  """ wurtzite lattice """
  from lada.crystal import Structure
  return Structure( 0.5, 0.5, 0,\
                    -0.866025, 0.866025, 0,\
                    0, 0, 1,\
                    scale=1, name='Wurtzite' )\
           .add_atom(0.5, 0.288675, 0, 'A')\
           .add_atom(0.5, -0.288675, 0.5, 'A')\
           .add_atom(0.5, 0.288675, 0.25, 'B')\
           .add_atom(0.5, -0.288675, 0.75, 'B')


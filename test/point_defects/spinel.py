
def lattice(x=0.387):
  from lada.crystal import Lattice

  lattice = Lattice()

  lattice.set_cell = (0, 0.5, 0.5),\
                     (0.5, 0, 0.5), \
                     (0.5, 0.5, 0) 

  # Octahedral sites
  lattice.add_sites = ( [5e0/8e0, 5e0/8e0, 5e0/8e0], "A" ),\
                      ( [5e0/8e0, 7e0/8e0, 7e0/8e0], "A" ),\
                      ( [7e0/8e0, 5e0/8e0, 7e0/8e0], "A" ),\
                      ( [7e0/8e0, 7e0/8e0, 5e0/8e0], "A" )
  # Tetrahedral sites
  lattice.add_sites = ( [   0,   0,   0], "B" ),\
                      ( [0.25,0.25,0.25], "B" )
  # Anion sites
  lattice.add_sites = ( [     x,      x,      x], "X" ),\
                      ( [     x,     -x,     -x], "X" ),\
                      ( [0.25-x, 0.25-x, 0.25-x], "X" ),\
                      ( [0.25-x, 0.25+x, 0.25+x], "X" ),\
                      ( [    -x,     -x,      x], "X" ),\
                      ( [    -x,      x,     -x], "X" ),\
                      ( [0.25+x, 0.25-x, 0.25+x], "X" ),\
                      ( [0.25+x, 0.25+x, 0.25-x], "X" )

  lattice.scale = 7.5

  lattice.find_space_group()
  assert len(lattice.space_group) != 0

  return lattice


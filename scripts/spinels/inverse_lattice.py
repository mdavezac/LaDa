
def lattice():
  from numpy import matrix, array
  from pylada import crystal

  lattice = crystal.Lattice()

  lattice.cell = matrix( [ [ 0, 0.5, 0.5 ], \
                           [ 0.5, 0, 0.5 ], \
                           [ 0.5, 0.5, 0 ] ] )

  # Manganese - Tetrahedral
  lattice.sites.append( crystal.Site(array([0, 0, 0], dtype="float64"), ["Mg"]) )
  lattice.sites.append( lattice.sites[0] )
  lattice.sites[0].pos = array( [   0,    0,    0], dtype="float64" )
  lattice.sites[1].pos = array( [0.25, 0.25, 0.25], dtype="float64" )

  # Aluminum - Octahedral
  lattice.sites.append( crystal.Site(array([5.0/8.0, 5.0/8.0, 5.0/8.0], dtype="float64"),\
                        ["Al", "Mg"]) )
  lattice.sites.extend( [ lattice.sites[2] for u in range(0,3) ] )
  lattice.sites[2].pos = array( [5.0/8.0, 5.0/8.0, 5.0/8.0], dtype="float64")
  lattice.sites[3].pos = array( [5.0/8.0, 7.0/8.0, 7.0/8.0], dtype="float64")
  lattice.sites[4].pos = array( [7.0/8.0, 5.0/8.0, 7.0/8.0], dtype="float64")
  lattice.sites[5].pos = array( [7.0/8.0, 7.0/8.0, 5.0/8.0], dtype="float64")

  # Oxygens
  x = 0.387
  lattice.sites.append( crystal.Site( array([x, x, x], dtype="float64"), ["O"]) )
  lattice.sites.extend( [ lattice.sites[6] for u in range(0,7) ] )
  lattice.sites[ 6].pos = array( [     x,      x,      x], dtype="float64" )
  lattice.sites[ 7].pos = array( [     x,     -x,     -x], dtype="float64" )
  lattice.sites[ 8].pos = array( [0.25-x, 0.25-x, 0.25-x], dtype="float64" )
  lattice.sites[ 9].pos = array( [0.25-x, 0.25+x, 0.25+x], dtype="float64" )
  lattice.sites[10].pos = array( [    -x,     -x,      x], dtype="float64" )
  lattice.sites[11].pos = array( [    -x,      x,     -x], dtype="float64" )
  lattice.sites[12].pos = array( [0.25+x, 0.25-x, 0.25+x], dtype="float64" )
  lattice.sites[13].pos = array( [0.25+x, 0.25+x, 0.25-x], dtype="float64" )

# trans = array([3.0/8.0, 3.0/8.0, 3.0/8.0], dtype="float64")
# for site in lattice.sites: site.pos -= trans

  lattice.scale = 7.5 # in bohr?
  lattice.find_space_group()

  return lattice

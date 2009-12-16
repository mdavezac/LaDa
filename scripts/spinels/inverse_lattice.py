
def lattice():
  from lada import crystal, atat

  lattice = crystal.Lattice()

  lattice.cell = atat.rMatrix3d( [ [ 0, 0.5, 0.5 ], \
                                   [ 0.5, 0, 0.5 ], \
                                   [ 0.5, 0.5, 0 ] ] )

  # Manganese
  lattice.sites.append( crystal.Site( (0, 0, 0) ) )
  lattice.sites[0].type = crystal.StringVector( [ "Mg" ] );
  lattice.sites.append( lattice.sites[0] )
  lattice.sites[0].pos = atat.rVector3d( (   0,    0,    0) )
  lattice.sites[1].pos = atat.rVector3d( (0.25, 0.25, 0.25) )

  # Aluminum
  lattice.sites.append( crystal.Site( (5.0/8.0, 5.0/8.0, 5.0/8.0) ) )
  lattice.sites[2].type = crystal.StringVector( ["Al", "Mg"] )
  lattice.sites.extend( [ lattice.sites[2] for u in range(0,3) ] )
  lattice.sites[2].pos = atat.rVector3d( [5.0/8.0, 5.0/8.0, 5.0/8.0] )
  lattice.sites[3].pos = atat.rVector3d( [5.0/8.0, 7.0/8.0, 7.0/8.0] )
  lattice.sites[4].pos = atat.rVector3d( [7.0/8.0, 5.0/8.0, 7.0/8.0] )
  lattice.sites[5].pos = atat.rVector3d( [7.0/8.0, 7.0/8.0, 5.0/8.0] )

  # Oxygens
  x = 0.387
  lattice.sites.append( crystal.Site( (x, x, x) ) )
  lattice.sites[6].type = crystal.StringVector( ["O"] )
  lattice.sites.extend( [ lattice.sites[6] for u in range(0,7) ] )
  lattice.sites[ 6].pos = atat.rVector3d( [     x,      x,      x] )
  lattice.sites[ 7].pos = atat.rVector3d( [     x,     -x,     -x] )
  lattice.sites[ 8].pos = atat.rVector3d( [0.25-x, 0.25-x, 0.25-x] )
  lattice.sites[ 9].pos = atat.rVector3d( [0.25-x, 0.25+x, 0.25+x] )
  lattice.sites[10].pos = atat.rVector3d( [    -x,     -x,      x] )
  lattice.sites[11].pos = atat.rVector3d( [    -x,      x,     -x] )
  lattice.sites[12].pos = atat.rVector3d( [0.25+x, 0.25-x, 0.25+x] )
  lattice.sites[13].pos = atat.rVector3d( [0.25+x, 0.25+x, 0.25-x] )

# for site in lattice.sites: 
#   site.pos -= atat.rVector3d( (3.0/8.0, 3.0/8.0, 3.0/8.0) )

  lattice.scale = 7.5

  lattice.find_space_group()
  assert (lattice.space_group) != 0

  return lattice

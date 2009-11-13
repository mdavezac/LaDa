#
#  Version: $Id$
#

class Mesh:

  def __init__( self, _start, _end, steps = 10, bandwidth = 2 ):
    from lada import crystal
    
    self.steps = steps
    self.bandwidth = bandwidth
    self.start = crystal.sStructure(_start)
    self.end = crystal.sStructure(_end)
    crystal.to_fractional( self.start )
    crystal.to_fractional( self.end )

  def generate( self ):
    from lada import crystal, atat
    
    id = atat.rMatrix3d( [[1,0,0], [0,1,0], [0,0,1]] )
    deformation = self.get_deformation()
    for n in range(0, self.steps+1):
      for a in range(n-self.bandwidth, n+self.bandwidth+1 ):

        structure = crystal.sStructure( self.start )
        structure.cell =   ( id + float(n) * deformation.cell ) \
                         * self.start.cell
        for i, atom in enumerate(structure.atoms):
          atom.pos += deformation.atoms[i].pos * float(a) 

        structure.name += "_%i_%i" % (n,a)

        crystal.to_cartesian(structure)
        yield structure

  def get_deformation( self ):
    from lada import crystal, atat

    id = atat.rMatrix3d( [[1,0,0], [0,1,0], [0,0,1]] )
    result = crystal.sStructure(self.end) 
    result.cell = ( self.end.cell * atat.inverse(self.start.cell) - id ) / float(self.steps)
    for i, atom in enumerate(result.atoms):
      atom.pos = ( self.end.atoms[i].pos - self.start.atoms[i].pos ) / float(self.steps)


    return result;

def main():
  from lada import crystal

  def create_lattice():

    from lada import crystal, atat

    lattice = crystal.Lattice()

    lattice.cell = atat.rMatrix3d( [ [ -1,  1,  1 ], \
                                     [  1, -1,  1 ], \
                                     [  1,  1, -1 ] ] )
    lattice.scale = 4.42

    lattice.sites.append( crystal.Site( (0, 0, 0) ) )
    lattice.sites[0].type = crystal.StringVector( [ "K", "Rb" ] );

    return lattice

  lattice = create_lattice()
  lattice.set_as_crystal_lattice()
  species = ("K", "Rb")
  start = crystal.read_poscar( species, "START" )
  end = crystal.read_poscar( species, "END" )

  for structure in Mesh(start, end, 4, 4).generate():
    print structure

if __name__ == "__main__":
  main()

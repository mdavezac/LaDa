class Mesh:

  def __init__( self, _start, _end, steps = 10, bandwidth = 2 ):
    from lada import crystal
    
    self.steps = steps
    self.bandwidth = bandwidth
    self.start = crystal.Structure(_start)
    self.end = crystal.Structure(_end)
    crystal.to_fractional( self.start )
    crystal.to_fractional( self.end )

  def generate( self ):
    from numpy import array as np_array, dot as np_dot
    from lada import crystal
    
    id = np_array( [[1,0,0], [0,1,0], [0,0,1]], dtype="float64" )
    deformation = self.get_deformation()
    for n in range(0, self.steps+1):
      for a in range(n-self.bandwidth, n+self.bandwidth+1 ):

        structure = crystal.Structure( self.start )
        structure.cell = np_dot(id + float(n) * deformation.cell, self.start.cell)
        for i, atom in enumerate(structure.atoms):
          atom.pos += deformation.atoms[i].pos * float(a) 

        structure.name += "_%i_%i" % (n,a)

        crystal.to_cartesian(structure)
        yield structure

  def get_deformation( self ):
    from numpy import array as np_array, dot as np_dot
    from lada import crystal

    id = np_array( [[1,0,0], [0,1,0], [0,0,1]], dtype="float64" )
    result = crystal.Structure(self.end) 
    result.cell = ( np_dot(self.end.cell, self.start.cell.I) - id ) / float(self.steps)
    for i, atom in enumerate(result.atoms):
      atom.pos = ( self.end.atoms[i].pos - self.start.atoms[i].pos ) / float(self.steps)


    return result;

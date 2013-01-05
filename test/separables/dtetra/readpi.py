def importPiStructure(file):
    """ load a Structure from Pi file

    We map the 1,2 indicator for atom type to 0,1
    We map the integer and halfinteger designation for BCC sites
       to even and odd integers

    filename and path are separated so that filename can be the
    structure name

    """
    import re
    from EquivStructure import Structure
    from numpy import array
    from scipy.linalg import inv

    # finds first instance of "NO."
    reNdec = re.compile('^\sNO\.')
    line = ""
    result = None
    for line in file:
      result = reNdec.match( line )
      if result is None: continue
      break
    if result is None: raise IOError 
    splitted = line.strip().split()
    # creates atomic species entries in structure
    N = eval(splitted[3])

    decoration = int( round( eval( splitted[4] ) ) )
    entries = [ int( (decoration >> n) %2 ) for n in range(N) ]
    cell = array( [ [ int( round( eval(x) ) ) for x in splitted[ 6: 9] ],
                    [ int( round( eval(x) ) ) for x in splitted[ 9:12] ],
                    [ int( round( eval(x) ) ) for x in splitted[12:15] ] ] )
    reciprocal=inv(cell)
    index = splitted[1]
    # now gets atomic positions
    reNdec = re.compile('^\sBASIS')
    current_pos = file.tell()
    for line in file:
      if reNdec.match( line ) is not None: break

    vectors = []
    oldline = line;

    splitted = line.strip().split()
    splitted = splitted[1:]
    for i in range( len(splitted) / 3 ):
      vectors.append(( int( round( eval( splitted[  i*3] ) ) ),
                       int( round( eval( splitted[i*3+1] ) ) ),
                       int( round( eval( splitted[i*3+2] ) ) ) ))

    if len( vectors ) < N:
      for line in file:
        splitted = line.strip().split()
        for i in range( len(splitted) / 3 ):
          vectors.append(( int( round( eval( splitted[  i*3] ) ) ),
                           int( round( eval( splitted[i*3+1] ) ) ),
                           int( round( eval( splitted[i*3+2] ) ) ) ))
        if len( vectors ) >= N: break
#   

    return Structure(vectors, entries,cell, reciprocal,index, 0.0)

def SFtoCE_Structure(structure):
    """ Converts EquivStructure.Structure to an atom.Structure
    """
    import Pylada
    from EquivStructure import Structure
    result = Pylada.Structure()
    dummy = Pylada.rMatrix3d(); dummy.diag( Pylada.rVector3d( [0.5, 0.5, 0.5 ] ) )
    result.cell = Pylada.make_rMatrix3d( structure.period ) * dummy
    
    for i in range( len( structure.entries ) ):
      result.atoms.append( Pylada.Atom( [ 0.5 * structure.vectors[i][0],
                                        0.5 * structure.vectors[i][1],
                                        0.5 * structure.vectors[i][2],
                                        structure.entries[i]*2-1 ] ) )
    result.scale = 1.0 
    return result
  

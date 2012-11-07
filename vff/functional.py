class Functional(object): 
  def __init__(self): 
    super(Functional, self).__init__()

    self._parameters = {}
    """ Holds vff parameters. """

  def __getitem__(self, index):
    """ VFF parameter access. 
    
        VFF parameters are composed of bond and angle parameters.
        
        Bond parameters can be accessed as:

        >>> functional['Au', 'Pd']
        array([2.62332, 21.6739, -112.0, 150.0, 0, 0])

        The order of the two species does not matter:

        >>> functional['Au', 'Pd'] is functional['Pd', 'Au']
        True

        The first item of the return array is the bond-length. The rest are
        the bond stretching parameters (second to 6th order).

        Angle parameters can be accessed the same way:

        >>> functional['Au', 'Pd', 'Cu']
        (-0.3333, -4.099, 9.3703)

        The first argument is the cosine of the equilibrium angle.
        The order of the first and last specie does not matter. The center of
        the angle is the central (second) specie.
    """
    if isinstance(index, str): index = index.split('-')
    index = '-'.join(sorted(str(u) for u in index))
    return self._parameters[index]

  def __setitem__(self, index, value):
    """ Adds/sets VFF parameters.
    
        VFF parameters are composed of bond and angle parameters.
        
        Bond parameters can be added/set as:

        >>> functional['Au', 'Pd'] = 2.62332, 21.6739, -112.0, 150.0, 0, 0

        The order of the two species does not matter. At this point, units
        should be based on angstroms. This is not checked in any way.

        The first item of the return array is the bond-length. The rest are
        the bond stretching parameters (second to 6th order).

        Angle parameters can be added/set in the same way:

        >>> functional['Au', 'Pd', 'Cu'] = -0.3333, -4.099, 9.3703

        The first argument is the cosine of the equilibrium angle.
        In the case of a perfect tetrahedron, it can be set from the string
        "tet". Furthermore, it can be set from a quantity signed in degrees
        or radians (see quantity_).

        The rest of the parameters are the bond-stretching parameters from
        2nd to 6th order. Their units should be based on angstroms though
        this is not checked in anyway.
    """
    from numpy import cos, array
    from quantities import radian
    from ..error import ValueError

    if isinstance(index, str): index = index.split('-')
    index = '-'.join(sorted(str(u) for u in index))
    # special case where the angle is given as "tet"
    maxsize = 7 if index.count('-') == 2 else 6
    if maxsize == 7 and isinstance(value[0], str):
      if value[0][:3].lower() != 'tet':
        raise ValueError( 'If a string, the first argument to angle '          \
                          'parameters should be "tet". ')
      value = [-1e0/3e0] + [u for u in value[1:]]
    # special case of a signed quantity.
    elif hasattr(value[0], 'rescale'):
      value = [cos(value[0].rescale(radian).magnitude)]                        \
              + [u for u in value[1:]]
    value = array(value).flatten()
    if len(value) < 2 or len(value) > maxsize:
      raise ValueError( 'Expects no less than two and no more than 6 '         \
                        'parameters.')
    self._parameters[index] = array( value.tolist()
                                     + [0]*(maxsize - len(value)) )

  def __call__(self, structure):
    """ Evaluates energy and forces on a structure. """
    from numpy import zeros
    from quantities import eV
    from . import build_tree
    # creates result structure.
    result = structure.copy()
    for atom in result: atom.gradient = zeros(3, dtype='float64')
    result.stress = zeros((3,3), dtype='float64')
    result.energy = 0e0
   
    scale2 = result.scale ** 2
    # creates tree and loop over structure.
    tree = build_tree(result)
    for node in tree: 
      result.energy += self._evaluate_bonds(node, scale2, result.cell,
                                              result.stress )
      result.energy += self._evaluate_angles( node, scale2, result.cell,
                                              result.stress )
    result.energy = result.energy / 16.0217733 * eV
    return result

  def _evaluate_bonds(self, node, scale2, cell, stress=None, strain=None,
                      K0=None):
    """ Evaluates bond-stretching modes. """
    from numpy import dot, sum, array, sqrt, outer
    bondparams = array([1e0, 5e-1/sqrt(3e0), 1e0/16e0, sqrt(3e0)/8e0, 1e0/640])
    gbondparams = array([ 1.5e0, 3e0*sqrt(3e0)/8.0, 3e0/16e0,
                          3e0*sqrt(3e0)/128e0, 0.00703125 ])

    energy = 0
    for endpoint, vector in node.sc_bond_iter():
      vector = dot(cell, vector) + endpoint.pos - node.pos
      params = self[node.type, endpoint.type]
      bond_length = params[0]

      e0 = sum(vector*vector) * scale2 / bond_length - bond_length
      mult = params[1:] * bondparams * e0
      energy += e0 * (mult[0] + e0 * ( mult[1]                                 \
                 + e0 * ( mult[2] + (mult[3] + e0 * mult[4]) ) ))

      if stress is not None: 
        mult = params[1:] * gbondparams * e0
        e0grad = 2e0 * scale2 / bond_length                                    \
                 * (mult[0] + e0 * ( mult[1]                                   \
                   + e0 * ( mult[2] + (mult[3] + e0 * mult[4]) ) ))
        hold = e0grad * (dot(strain, vector) if strain is not None else vector)
        node.center.gradient -= hold
        endpoint.center.gradient += hold
  
        matrix = outer(vector, vector)
        stress += e0grad * 0.5 * (dot(matrix, K0) if K0 is not None else matrix)
                
    return energy * 3e0 / 8e0 

  def _evaluate_angles( self, node, scale2, cell, stress=None, strain=None,
                        K0=None ):
    """ Evaluates bond-angle and bond-bending modes. """
    from numpy import dot, sum, array, sqrt, outer
    bondparams = array([1e0, 5e-1/sqrt(3e0), 1e0/16e0, sqrt(3e0)/8e0, 1e0/640])
    gbondparams = array([ 0.75e0, 3e0*sqrt(3e0)/16.0, 3e0/32e0,
                          3e0*sqrt(3e0)/256e0, 0.003515625 ])

    energy = 0
    for (A,dA), (B, dB) in node.angle_iter():
      vA = dot(cell, dA) + A.pos - node.pos
      vB = dot(cell, dB) + B.pos - node.pos
      paramsA = self[node.type, A.type]
      paramsB = self[node.type, B.type]
      paramsAB = self[A.type, node.type, B.type]
      lengthA, lengthB = paramsA[0], paramsB[0]
      gamma, sigma =  paramsAB[0], paramsAB[1]
      mean_length = sqrt(lengthA * lengthB)

      e0 = sum(vA*vA) * scale2 / lengthA - lengthA                             \
            + sum(vB*vB) * scale2 / lengthB - lengthB
      e1 = dot(vA, vB) * scale2 / mean_length - mean_length * gamma 

      # bond-bending
      mult = paramsAB[2:] * bondparams * e1
      energy += e1 * (mult[0] + e1 * ( mult[1]                                 \
                 + e1 * ( mult[2] + (mult[3] + e1 * mult[4]) ) ))

      if stress is not None: 
        mult = paramsAB[2:] * gbondparams * e1
        e1grad = 2e0 * scale2 / mean_length                                    \
                 * (mult[0] + e1 * ( mult[1]                                   \
                   + e1 * ( mult[2] + (mult[3] + e1 * mult[4]) ) ))
        hold0 = e1grad * (dot(strain, vA) if strain is not None else vA)
        hold1 = e1grad * (dot(strain, vB) if strain is not None else vB)
        node.center.gradient -= hold0 + hold1
        A.center.gradient += hold1
        B.center.gradient += hold0
  
        matrix = outer(vA, vB)
        matrix += matrix.T
        if K0 is not None: matrix = dot(matrix, K0)
        stress += e1grad * 0.5 * matrix

      # bond angle 
      energy += e0 * e1 * sigma
      
      if stress is not None: 
        defA = dot(strain, vA) if strain is not None else vA
        defB = dot(strain, vB) if strain is not None else vB
        hold0 = 1.5 * e1 * sigma / lengthA * scale2 * defA                     \
                + 0.75 * e0 * sigma / mean_length * scale2 * defB 
        hold1 = 1.5 * e1 * sigma / lengthB * scale2 * defB                     \
                + 0.75 * e0 * sigma / mean_length * scale2 * defA


        node.center.gradient -= hold0 + hold1
        A.center.gradient += hold0
        B.center.gradient += hold1

        matrix = outer(vA, vB)
        matrix = 2e0 * e1 * (outer(vA, vA)/lengthA + outer(vB, vB)/lengthB)    \
                 + e0 / mean_length * (matrix + matrix.T)
        if K0 is not None: matrix = dot(matrix, K0)
        stress += matrix * 0.375 * sigma * scale2

    return energy * 3e0 / 8e0


class Functional(object): 
  def __init__(self): 
    super(Functional, self).__init__()

    self._parameters
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
    if index.count('-') == 2 and isinstance(value[0], str):
      if value[0][:3].lower() != 'tet':
        raise ValueError( 'If a string, the first argument to angle '          \
                          'parameters should be "tet". ')
        value = [-1e0/3e0] + [u for u in value[1:]]
    # special case of a signed quantity.
    elif hasattr(value[0], 'rescale'):
      value = [cos(value[0].rescale(radian).magnitude)]                        \
              + [u for u in value[1:]]
    value = array(value).flatten()
    if len(value) < 2 or len(value) > 6:
      raise ValueError( 'Expects no less than two and no more than 6 '         \
                        'parameters.')
    self._parameters[index] = array(value.tolist() + [0]*(6 - len(value)))

# def __call__(self, structure):
#   """ Evaluates energy and forces on a structure. """
#   from numpy import zeros
#   from . import build_tree
#   # creates result structure.
#   result = structure.copy()
#   for atom in result: result.gradient = zeros(3, dtype='float64')
#   result.stress = zeros((3,3), dtype='float64')
#   result.energy = 0e0
#  
#   scale2 = result.scale ** 2
#   # creates tree and loop over structure.
#   tree = build_tree(result)
#   for node in tree: 
#     energy, dself._evaluate_bonds(node, scale)
#     length = 
#     for bond in node:

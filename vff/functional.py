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
    from collections import Sequence
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
    from numpy import cos
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
#   # creates result structure.
#   result = structure.copy()
#   for atom in result: result.gradient = zeros(3, dtype='float64')
#   result.stress = zeros((3,3), dtype='float64')
#   result.energy = 0e0
#  
#   scale2 = result.scale ** 2
#   # creates tree and loop over structure.
#   tree = self.build_tree(result)
#   for node in tree: 
#     energy, dself._evaluate_bonds(node, scale)
#     length = 
#     for bond in node:

      

  @staticmethod
  def build_tree(structure, overlap=1.2, **kwargs):
    """ Build first neighbor tree. """
    from numpy import array, argsort, dot, sum
    from numpy.linalg import inv
    from quantities import angstrom
    from ..crystal import periodic_dnc
    from ..error import ValueError
    from . import Node

    if hasattr(overlap, 'rescale'):
      overlap = float(overlap.rescale(angstrom)) / structure.scale

    # net: map of atom id to the nodes in the the first neighbor net which
    #      wraps them.
    net = {}
    for atom in structure: net[id(atom)] = Node(atom)
    
    # invcell: transform from cartesian to fractional coordinates.
    invcell = inv(structure.cell)
    # boxes: list of divide and conquer boxes.
    boxes = periodic_dnc(structure, overlap=overlap, **kwargs)
    # loop over boxes.
    for box in boxes:
      if len(box) < 5: raise ValueError('Box does not contain enough atoms.')
      # position: All positions in box as a numpy array, for faster treatment.
      positions = array([a.pos + b for a, b, c in box])
      # Loop over each node in the net in that particular box.
      for position, (atom, vector, inbox) in zip(positions, box):
        # If atom is outside box, loop.
        if not inbox: continue 
        # node: Node which wraps current atom
        node = net[id(atom)]
        if len(node) == 4: continue # already found all bonds.

        # distances: distance from this atom to all others in box.
        distances = sum((positions - position)**2, axis=1)
        # sortees: indices of four nearest neighbors in box. 
        sortees = argsort(distances)[1:5]
        
        # Now adds sortee as bonds.
        for i in sortees: 
          frac = dot(invcell, box[i][1] - vector)
          othernode = net[id(box[i][0])]
          if len(othernode) == 4: continue
          node.link(othernode, frac)
          if len(node) == 4: break
    # returns net as a list.
    return net.values()

from lada import is_interactive
if is_interactive:
  def plot_tree(structure, overlap=0.5):
    """ Plots a tree using mayavi. """
    from enthought.mayavi import mlab
    from numpy import array, sum, any, dot, concatenate, argmin
    from lada.vff.functional import Functional
    tree = Functional.build_tree(structure, overlap=overlap)
 
    # find original positions
    ids, originals = [], []
    for i, node in enumerate(tree):
      if id(node.center) not in ids: 
        ids.append(id(node.center))
        originals.append(node.pos)
    originals = array(originals)
    # now adds neighbors, including those outside the unit cell
    positions = originals.copy()
    others, connections = [], []
    for i, node in enumerate(tree): 
      for neighbor, vector in node:
        position = neighbor.pos + dot(structure.cell, vector)
        dist = sum((positions-position)**2, axis=1) 
        if any(dist < 1e-8): connections.append( (i, argmin(dist)) )
        else: 
          positions = concatenate((positions, position[None, :]))
          connections.append((i, positions.shape[0]-1))
           
        if any(sum((position[None,:]-originals)**2, axis=1) < 1e-8): continue
        if len(others)                                                           \
           and any(sum((position[None,:]-others)**2, axis=1) < 1e-8): continue
        others.append(position)
    others = array(others)
 
    a0, a1, a2 = structure.cell.T
    box = [ [0, 0, 0], a0, a0+a1, a1, [0,0,0], a2, a2+a0, a2+a1+a0, a2+a1, a2,  \
            a2 + a0, a0, a0+a1, a0+a1+a2, a1+a2, a1 ] 
    mlab.plot3d(*array(box).T, tube_radius=0.01, color=(0.8, 0.8, 0.8))
 
    mlab.points3d(*originals.T, color=(0.8, 0.5, 0.5), scale_factor=0.2)
    mlab.points3d(*others.T, color=(0.5, 0.8, 0.5), scale_factor=0.1)
    # add connections
    pts = mlab.points3d(*positions.T, scale_factor=0)
    pts.mlab_source.dataset.lines = array(connections)
    tube = mlab.pipeline.tube(pts, tube_radius=0.01)
    tube.filter.radius_factor = 1.
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    mlab.pipeline.surface(tube, color=(0.3, 0.3, 0.3))

  def plot_box(cell, tube_radius=0.01, color=None):
    """ Adds rectangular box to mayavi. """
    from enthought.mayavi import mlab 
    from numpy import array
    a0, a1, a2 = cell.T
    box = array([[0, 0, 0], a0, a0+a1, a1, a2, a2+a0, a2+a1+a0, a2+a1])
    connections = array([ [0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [4, 5], 
                          [5, 6], [6, 7], [7, 4], [1, 5], [7, 3], [6, 2] ])

    pts = mlab.points3d(*box.T, scale_factor=0)
    pts.mlab_source.dataset.lines = connections
  
    # Use a tube fiter to plot tubes on the link, varying the radius with the
    # scalar value
    tube = mlab.pipeline.tube(pts, tube_radius=tube_radius)
    tube.filter.radius_factor = 1.
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    if color is None: color = (0.8, 0.8, 0)
    mlab.pipeline.surface(tube, color=color)





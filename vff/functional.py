class Functional(object):
  def __init__(self): 
    super(Functional, self).__init__()

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
          frac = dot(invcell, positions[i] - position)
          othernode = net[id(box[i][0])]
          node.link(othernode, frac)
          if len(node) == 4: break
    # returns net as a list.
    return net.values()

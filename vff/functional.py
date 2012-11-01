class Functional(object): pass
# def __init__(self): 
#   super(Functional, self).__init__()

# @staticmethod
# def build_tree(structure, overlap=1.2, **kwargs):
#   """ Build first neighbor tree. """
#   from numpy import array, argsort, dot, sum
#   from numpy.linalg import inv
#   from quantities import angstrom
#   from ..crystal import periodic_dnc
#   from ..error import ValueError
#   from . import Node

#   if hasattr(overlap, 'rescale'):
#     overlap = float(overlap.rescale(angstrom)) / structure.scale

#   # net: map of atom id to the nodes in the the first neighbor net which
#   #      wraps them.
#   net = {}
#   for atom in structure: net[id(atom)] = Node(atom)
#   
#   # invcell: transform from cartesian to fractional coordinates.
#   invcell = inv(structure.cell)
#   # boxes: list of divide and conquer boxes.
#   boxes = periodic_dnc(structure, overlap=overlap, **kwargs)
#   # loop over boxes.
#   for box in boxes:
#     if len(box) < 5: raise ValueError('Box does not contain enough atoms.')
#     # position: All positions in box as a numpy array, for faster treatment.
#     positions = array([a.pos + b for a, b, c in box])
#     # Loop over each node in the net in that particular box.
#     for position, (atom, vector, inbox) in zip(positions, box):
#       # If atom is outside box, loop.
#       if not inbox: continue 
#       # node: Node which wraps current atom
#       node = net[id(atom)]
#       if len(node) == 4: continue # already found all bonds.

#       # distances: distance from this atom to all others in box.
#       distances = sum((positions - position)**2, axis=1)
#       # sortees: indices of four nearest neighbors in box. 
#       sortees = argsort(distances)[1:5]
#       
#       # Now adds sortee as bonds.
#       for i in sortees: 
#         frac = dot(invcell, box[i][1] - vector)
#         othernode = net[id(box[i][0])]
#         if len(othernode) == 4: continue
#         print frac
#         node.link(othernode, frac)
#         if len(node) == 4: break
#   # returns net as a list.
#   return net.values()


from lada import is_interactive
if is_interactive:
  def plot_tree(structure):
    """ Plots a tree using mayavi. """
    pass
#   from mayavi import mlab
#   from numpy import array, sum, any, dot
#   from lada.vff.functional import Functional
#   tree = Functional.build_tree(structure, overlap=0.5)
#
#   # find original positions
#   ids, originals = [], []
#   for i, node in enumerate(tree):
#     if id(node.center) not in ids: 
#       ids.append(id(node.center))
#       originals.append(node.pos)
#   originals = array(originals)
#   # now adds neighbors, including those outside the unit cell
#   others, lines = [], []
#   for node in tree: 
#     for neighbor, vector in node:
#       position = neighbor.pos + dot(structure.cell, vector)
#       lines.append([node.pos, position])
#       if any(sum((position[None,:]-originals)**2, axis=1) < 1e-8): continue
#       if len(others)                                                           \
#          and any(sum((position[None,:]-others)**2, axis=1) < 1e-8): continue
#       others.append(position)
#   others = array(others)
#
#   a0, a1, a2 = structure.cell.T
#   box = [ [0, 0, 0], a0, a0+a1, a1, [0,0,0], a2, a2+a0, a2+a1+a0, a2+a1, a2,  \
#           a2 + a0, a0, a0+a1, a0+a1+a2, a1+a2, a1 ] 
#   mlab.plot3d(*array(box).T, tube_radius=0.01, color=(0.8, 0.8, 0.8))
#
#   for a, b in lines:
#     l = array([a, b])
#     mlab.plot3d(*l.T, color=(0.8, 0.5, 0.5), tube_radius=0.01)
#
#   mlab.points3d(*originals.T, color=(0.8, 0.5, 0.5), scale_factor=0.2)
#   mlab.points3d(*others.T, color=(0.5, 0.8, 0.5), scale_factor=0.1)

  def plot_box(cell, tube_radius=0.01, color=None):
    """ Adds rectangular box to mayavi. """
    pass
#   from mayavi import mlab 
#   from numpy import array
#   a0, a1, a2 = cell.T
#   box = array([[0, 0, 0], a0, a0+a1, a1, a2, a2+a0, a2+a1+a0, a2+a1])
#   connections = array([ [0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [4, 5], 
#                         [5, 6], [6, 7], [7, 4], [1, 5], [7, 3], [6, 2] ])

#   pts = mlab.points3d(*box.T, scale_factor=0)
#   pts.mlab_source.dataset.lines = connections
# 
#   # Use a tube fiter to plot tubes on the link, varying the radius with the
#   # scalar value
#   tube = mlab.pipeline.tube(pts, tube_radius=tube_radius)
#   tube.filter.radius_factor = 1.
#   tube.filter.vary_radius = 'vary_radius_by_scalar'
#   if color is None: color = (0.8, 0.8, 0)
#   mlab.pipeline.surface(tube, color=color)





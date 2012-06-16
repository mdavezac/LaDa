class Spin(object):
  __slots__ = ['position', 'flavor', 'sublattice']
  def __init__(self, position=None, flavor=0, sublattice=0):
    """ Sets up a spin, """
    self.position 
class Cluster(object):
  __slots__ = ['interactions', 'flavor', 'site']
  def __init__(self, interactions=None, flavor=0, site=-1):
    from copy import deepcopy
    super(Cluster, self).__init__()
    self.interactions = deepcopy(interactions)
    """ Interactions involved in cluster. """
    self.flavor = flavor
    """ Flavor 

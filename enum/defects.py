class Iterator(object):

  def __init__(self, template, *args): 
    super(Iterator, self).__init__()
    self.template = template.copy()
    self.others = args
    self.reset()
    
  def __iter__(self): return self

  def reset(self):
    from numpy import count_nonzero, ones, logical_and
    from ..error import ValueError
    from .cppwrappers import FCIterator

    # creates list of color iterators, as per Hart, Nelson, Forcade.
    self.iterators = []
    current_allowed = ones(len(self.template), dtype='bool')
    for n, color, mask in self.others:
      allowed = logical_and(current_allowed, mask)
      length = count_nonzero(allowed)
      if length < n:
        raise ValueError( 'Could not create iterator, '                        \
                          'concentrations do not add up.')
      self.iterators.append(FCIterator(length, n))
      current_allowed[allowed] = [False] * n + [True] * (length-n)
    # necessary cos we bypass calling next for the first time in this
    # instance's next.
    for iter in self.iterators[1:]: iter.next()

    self.x = self.template.copy()

  def next(self):
    from numpy import ones, logical_and
    from ..error import internal
      
    self.x[:] = self.template
    mask = ones(len(self.x), dtype='bool')
    donext = True
    for iter, (n, color, cmask) in zip(self.iterators, self.others):
      if donext: 
        try: bitstring = iter.next()
        except StopIteration:
          iter.reset()
          try: bitstring = iter.next()
          except StopIteration: 
            raise internal('Cannot iterate over type {0}'.format(color))
        else: donext = False
      else: bitstring = iter.yielded
      change_color = logical_and(cmask, mask)
      change_color[change_color] = bitstring
      self.x[change_color] = color
      mask[change_color] = False
    if donext: raise StopIteration
    return self.x

# def defects(lattice, n, defects):
#   """ Generates defects on a lattice """
#   from ..error import ValueError
#   from .transforms import Transforms
#   from . import hf_groups

#   # sanity checks.
#   # Sites without defects should be spectating.
#   # No site with more than two types.
#   # One defect type per site.
#   # Sets defect type to second flavor.
#   lattice = lattice.copy()
#   if any(len(u.type) > 2 for u in lattice):
#     raise ValueError('Lattice should be a pseudo-binary.')
#   transforms = Transforms(lattice)
#   lattice = Transforms.lattice
#   for site in lattice:
#     if lattice.nbflavors == 1: continue
#     if len(set(defects.keys()).intersection(set(site.type))) > 2:
#       raise ValueError('More than one defect per site.')
#     if len(set(site.type) - set(defects.keys())) > 2:
#       raise ValueError('Multinary site without defects.')
#     for key in defects.keys():
#       if key in site.type and site.type == 2 and site.type[1] != key:
#         a, b = site.type[0], site.type[1]
#         site.type[0], site.type[1] = b, a
#     for key in defects.keys():
#       if key in site.type and site.type == 2 and site.type[1] != key:
#         raise ValueError('defects are on same lattice site.')
#   lattice = Transforms.lattice

#   # Find atom with the minimum number of inequivalent lattice sites.
#   # Its position will be fixed
#   fixedways = {}
#   for key in defects:
#     fixedways[key] = [ u for u in transforms.lattice                           \
#                        if u.nbflavors > 1 and key in u.type                    \
#                           and u.asymmetric ] 
#   fixedatom = min( (item for item in fixedways.iteritems()),                   \
#                    key = lambda x: len(x[1]) )

#   defects = defects.copy()
#   defects[fixedatom[0]] -= 1
#   if defects[fixedatom[0]] == 0: del defects[fixedatom[0]]

#   # loop over sizes
#   for hfgroups in hf_groups(lattice, n):
#     # loop over Hart-Forcade equivalent supercells
#     for hfgroup in hfgroups:
#       # actual results
#       ingroup = []
#       # stuff we do not want to see again
#       outgroup = set()
#   
#       # translation operators
#       translations = transforms.translations(hfgroup[0][0])
#   
#       # Creates argument to ndimensional iterator.
#       first = [u.index for u in fixedatom[1]]
#       args = []
#       size = hfgroup[0][0].size
#       for key in defects:
#         args += [size * len([u for u in lattice if key in u.type])]

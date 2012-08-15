class Iterator(object):

  def __init__(self): 
    from numpy import array
    super(Iterator, self).__init__()
    self.x = array()
  def next():
    result = 0


def defects(lattice, n, defects):
  """ Generates defects on a lattice """
  from ..error import ValueError
  from .transforms import Transforms

  # sanity checks.
  # Sites without defects should be spectating.
  # No site with more than two types.
  # One defect type per site.
  # Sets defect type to second flavor.
  lattice = lattice.copy()
  if any(len(u.type) > 2 for u in lattice):
    raise ValueError('Lattice should be a pseudo-binary.')
  transforms = Transforms(lattice)
  lattice = Transforms.lattice
  for site in lattice:
    if lattice.nbflavors == 1: continue
    if len(set(defects.keys()).intersection(set(site.type))) > 2:
      raise ValueError('More than one defect per site.')
    if len(set(site.type) - set(defects.keys())) > 2:
      raise ValueError('Multinary site without defects.')
    for key in defects.keys():
      if key in site.type and site.type == 2 and site.type[1] != key:
        a, b = site.type[0], site.type[1]
        site.type[0], site.type[1] = b, a
    for key in defects.keys():
      if key in site.type and site.type == 2 and site.type[1] != key:
        raise ValueError('defects are on same lattice site.')
  lattice = Transforms.lattice

  # Find atom with the minimum number of inequivalent lattice sites.
  # Its position will be fixed
  fixedways = {}
  for key in defects:
    fixedways[key] = [ u for u in transforms.lattice                           \
                       if u.nbflavors > 1 and key in u.type                    \
                          and u.asymmetric ] 
  fixedatom = min( (item for item in fixedways.iteritems()),                   \
                   key = lambda x: len(x[1]) )

  defects = defects.copy()
  defects[fixedatom[0]] -= 1
  if defects[fixedatom[0]] == 0: del defects[fixedatom[0]]

  # loop over sizes
  for hfgroups in hf_groups(lattice, n):
    # loop over Hart-Forcade equivalent supercells
    for hfgroup in hfgroups:
      # actual results
      ingroup = []
      # stuff we do not want to see again
      outgroup = set()
  
      # translation operators
      translations = transforms.translations(hfgroup[0][0])
  
      # Creates argument to ndimensional iterator.
      first = [u.index for u in fixedatom[1]]
      args = []
      size = hfgroup[0][0].size
      for key in defects:
        args += [size * len([u for u in lattice if key in u.type])]

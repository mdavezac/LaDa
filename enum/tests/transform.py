def get_cell(n=5):
  from numpy.random import randint
  from numpy.linalg import det
  cell = randint(2*n, size=(3,3)) - n
  while abs(det(cell)) < 1e-8:
    cell = randint(2*n, size=(3,3)) - n
  return cell

def test_translations():
  from numpy import abs, all
  from lada.crystal import binary, supercell, HFTransform
  from lada.enum import Transforms

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  for u in xrange(10):
    # create random structure
    structure = supercell(lattice, get_cell(5))
    hft = HFTransform(lattice, structure)
 
    # these are all the translations
    translations = Transforms(lattice).translations(hft)
    assert translations.shape == (len(structure)//len(lattice) - 1, len(structure))
    # compute each translation and gets decorations
    for atom in structure:
      if atom.site != 0: continue
      # create translation
      trans = atom.pos - lattice[0].pos
      if all(abs(trans) < 1e-8): continue
      # figure out its index
      index = hft.index(trans) - 1
      for site in structure:
        pos = site.pos-lattice[site.site].pos 
        i = hft.index(pos, site.site)
        j = hft.index(pos + trans, site.site)
        assert translations[index, i] == j

def test_firstisid():
  """ Assumption is made in Transforms.transformations """
  from numpy import abs, all, identity
  from lada.crystal import binary, supercell, space_group

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  for u in xrange(10):
    # create random structure
    structure = supercell(lattice, get_cell(5))
    while len(structure) > 1: structure.pop(-1)
    assert len(structure) == 1
    sg = space_group(structure)[0]
    assert all(abs(sg[:3]- identity(3, dtype='float64'))<1e-8)
    assert all(abs(sg[3])<1e-8)

def test_rotations():
  from numpy import all, dot, zeros
  from numpy.linalg import inv
  from lada.crystal import binary, supercell, HFTransform, space_group,        \
                           which_site
  from lada.enum import Transforms

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']
  sg = space_group(lattice)
  invcell = inv(lattice.cell)

  def get_cells(n):
    for i in xrange(1, n):
      yield [[i, 0, 0], [0, 0.5, 0.5], [0, -0.5, 0.5]]
    for i in xrange(1, n):
      yield [[i, 0, 0], [0, i, 0], [0, 0, 1]]
    for i in xrange(1, n):
      yield [[i, 0, 0], [0, i, 0], [0, 0, i]]
  for cell in get_cells(8):
    # create random structure
    structure = supercell(lattice, cell)
    hft = HFTransform(lattice, structure)
 
    # these are all the translations
    transforms = Transforms(lattice)
    permutations = transforms.transformations(hft)
    assert permutations.shape == (len(sg) - 1, len(structure))
    operations = transforms.invariant_ops(structure)
    assert any(operations) 

    # compute each translation and gets decorations
    for index, (op, isgood) in enumerate(zip(sg[1:], operations)):
      if not isgood: continue
      # Create rotation and figure out its index
      permutation = zeros(len(structure), dtype='int') - 1
      for atom in structure:
        pos = dot(op[:3], atom.pos) + op[3]
        newsite = which_site(pos, lattice, invcell)
        i = hft.index(atom.pos - lattice[atom.site].pos, atom.site)
        j = hft.index(pos - lattice[newsite].pos, newsite)
        permutation[i] = j
      assert all(permutation == permutations[index])

def test_lattice():
  """ Tests lattice enhancements. """
  from numpy import all
  from lada.crystal import binary, A2BX4
  from lada.enum import Transforms

 #lattice = binary.zinc_blende()
 #lattice[0].type = ['Si', 'Ge']
 #lattice[1].type = ['Si', 'Ge']
 #transforms = Transforms(lattice)
 #assert len([u for u in transforms.lattice if u.asymmetric]) == 1
 #assert transforms.lattice[0].asymmetric
 #assert transforms.lattice[0].equivto == 0
 #assert transforms.lattice[0].nbflavors == 2
 #assert transforms.lattice[0].index == 0
 #assert not transforms.lattice[1].asymmetric
 #assert transforms.lattice[1].equivto == 0
 #assert transforms.lattice[1].nbflavors == 2
 #assert transforms.lattice[1].index == 1
 #assert all(all(a == b) for a, b in zip(transforms.flavors, (range(1), range(1))))
 #assert all(not hasattr(atom, 'nbflavors') for atom in lattice)

 #lattice = binary.zinc_blende()
 #lattice[0].type = ['Si', 'Ge']
 #lattice[1].type = ['Si', 'Ge', 'C']
 #transforms = Transforms(lattice)
 #assert len([u for u in transforms.lattice if u.asymmetric]) == 2
 #assert transforms.lattice[0].asymmetric
 #assert transforms.lattice[0].equivto == 0
 #assert transforms.lattice[0].nbflavors == 2
 #assert transforms.lattice[0].index == 0
 #assert transforms.lattice[0].asymmetric
 #assert transforms.lattice[1].equivto == 1
 #assert transforms.lattice[1].nbflavors == 3
 #assert transforms.lattice[1].index == 1
 #assert all(all(a == b) for a, b in zip(transforms.flavors, (range(1), range(2))))

  lattice = A2BX4.b5() 
  for atom in lattice: 
    if atom.type in ['A', 'B']: atom.type = 'A', 'B'
  transforms = Transforms(lattice)
  print transforms.lattice
 #assert len([u for u in transforms.lattice if u.asymmetric]) == 3
 #assert all([transforms.lattice[i].asymmetric for i in [0, 4, 6]])
 #assert all([not transforms.lattice[i].asymmetric for i in range(1, 4) + [5] + range(7, 14)])
 #assert all([transforms.lattice[i].equivto == 0 for i in range(4)])
 #assert all([transforms.lattice[i].equivto == 4 for i in range(4, 6)])
 #assert all([transforms.lattice[i].equivto == 6 for i in range(6, 14)])
 #assert all([transforms.lattice[i].nbflavors == 2 for i in range(4)])
 #assert all([transforms.lattice[i].nbflavors == 2 for i in range(4, 6)])
 #assert all([transforms.lattice[i].nbflavors == 1 for i in range(6, 14)])
 #assert all([transforms.lattice[i].index == i for i in range(6)])
 #assert all([not hasattr(transforms.lattice[i], 'index') for i in range(6, 14)])

  lattice[0], lattice[-1] = lattice[-1], lattice[0]
  print lattice
  transforms = Transforms(lattice)
  print transforms.lattice

def test_toarray():
  """ Tests label exchange """
  from random import choice
  from numpy import all, dot, zeros
  from numpy.linalg import inv
  from lada.crystal import binary, supercell, HFTransform, space_group,        \
                           which_site
  from lada.enum import Transforms

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']
  sg = space_group(lattice)
  invcell = inv(lattice.cell)
  transforms = Transforms(lattice)
  lattice = transforms.lattice

  for u in xrange(11):
    structure = supercell(lattice, get_cell())
    for atom in structure: atom.type = choice(atom.type)
    hft = HFTransform(lattice, structure) 
    a = transforms.toarray(hft, structure)
    b = zeros(len(structure), dtype='int')
    for atom in structure:
      site = lattice[atom.site]
      b[hft.index(atom.pos-site.pos, atom.site)]                               \
        = site.type.index(atom.type) + 1
    assert all(a == b)


if __name__ == '__main__':
# test_translations()
# test_firstisid()
# test_rotations()
  test_lattice()
# test_toarray()


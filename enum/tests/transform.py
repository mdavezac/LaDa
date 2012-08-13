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
  from lada.crystal import binary
  from lada.enum import Transforms

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']
  transforms = Transforms(lattice)
  assert len([u for u in transforms.lattice if u.asymmetric]) == 1
  assert transforms.lattice[0].asymmetric
  assert transforms.lattice[0].equivto == 0
  assert not transforms.lattice[1].asymmetric
  assert transforms.lattice[1].equivto == 0
  assert all(all(a == b) for a, b in zip(transforms.flavors, (range(1), range(1))))

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']
  transforms = Transforms(lattice)
  assert len([u for u in transforms.lattice if u.asymmetric]) == 2
  assert transforms.lattice[0].asymmetric
  assert transforms.lattice[0].equivto == 0
  assert transforms.lattice[0].asymmetric
  assert transforms.lattice[1].equivto == 1
  assert all(all(a == b) for a, b in zip(transforms.flavors, (range(1), range(2))))

# def test_labelexchange():
#   """ Tests label exchange """
#   from random import choice
#   from numpy import all, dot, zeros
#   from numpy.linalg import inv
#   from lada.crystal import binary, supercell, HFTransform, space_group,        \
#                            which_site
#   from lada.enum import Transforms

#   lattice = binary.zinc_blende()
#   lattice[0].type = ['Si', 'Ge']
#   lattice[1].type = ['Si', 'Ge', 'C']
#   sg = space_group(lattice)
#   invcell = inv(lattice.cell)

#   for u in xrange(1):
#     structure = supercell(lattice, get_cell(u))
#     transforms = Transforms(lattice)
#     assert len([u for u in transforms.lattice if u.asymmetric]) == 2
#     for atom in structure: atom.type = choice(atom.type)
#     hft = HFTransform(lattice) 

if __name__ == '__main__':
# test_translations()
# test_firstisid()
# test_rotations()
  test_lattice()

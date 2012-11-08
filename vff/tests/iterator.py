def test_single_counting():
  from lada.crystal import binary, supercell
  from lada.vff import build_tree
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[4, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = build_tree(a, overlap=0.5)
  
  n = 0
  for center in b:
    for endpoint, vector in center.sc_bond_iter():
      n += 1
      for other, v in endpoint.sc_bond_iter(): assert other is not center
      assert id(center) in [id(c) for c, v in endpoint]
  assert n == 2 * len(a)

def test_angle():
  from lada.crystal import binary, supercell
  from lada.vff import build_tree

  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[4, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = build_tree(a, overlap=0.5)

  for center in b:
    ids = [id(u.center) for u, d in center]
    angles = set([(id(u.center), id(v.center)) for (u, d), (v, d) in center.angle_iter()])
    for i, ida in enumerate(ids[:-1]):
      for idb in ids[i+1:]:
        if (ida, idb) in angles: assert (idb, ida) not in angles
        else: assert (idb, ida) in angles
    assert len(angles) == 6

if __name__ == '__main__':
  test_single_counting()
  test_angle()

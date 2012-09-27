def test():
  from numpy import all, abs, dot, identity
  from lada.dftcrystal import Crystal, External
  from lada.crystal import which_site

  orig = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  b = External(copy=orig.eval())
  assert all(abs(orig.eval().cell - b.initial.cell) < 1e-8)
  assert len(b.initial) == 2
  assert repr(b) == repr( eval(repr(b), {'External': External}) )
  assert len(b.print_input().splitlines()) == 2
  assert b.print_input().splitlines()[0] == 'EXTERNAL'
  assert b.print_input().splitlines()[1] == 'END EXTERNAL'

  # check that symmetries as extracted are ok.
  for atom in b.initial:
    for op in orig.symmetry_operators:
      assert all(abs(identity(3, dtype='float64') - dot(op[:3].T, op[:3])) < 1e-8)
      assert which_site(dot(op[:3], atom.pos) + op[3], b.initial) != -1

  # Check that running structure through crystal works
  c = b.eval()
  assert all(abs(c.cell - b.cell) < 1e-8)
  assert len(c) == 2
  assert all(abs(abs(c[0].pos) - [0.67875]*3) < 1e-8)
  assert all(abs(abs(c[1].pos) - [0.67875]*3) < 1e-8)
  assert all(abs(c[0].pos + c[1].pos) < 1e-8)

if __name__ == '__main__':
  test()

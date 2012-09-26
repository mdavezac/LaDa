def test():
  from numpy import all, abs
  from lada.dftcrystal import Crystal, External
  from lada.crystal import write

  orig = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  b = External(copy=orig.eval())
  assert all(abs(orig.eval().cell - b.initial.cell) < 1e-8)
  assert len(b.initial) == 2
  assert repr(b) == repr( eval(repr(b), {'External': External}) )
  assert len(b.print_input().splitlines()) == 2
  assert b.print_input().splitlines()[0] == 'EXTERNAL'
  assert b.print_input().splitlines()[1] == 'END EXTERNAL'
  print write.crystal(b.initial, None)
  print b.eval()

if __name__ == '__main__':
  test()

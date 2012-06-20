def test_crystal():
  from numpy import all, abs
  from lada.dftcrystal.geometry import Crystal
  a = Crystal()
  a.raw = '0 0 0\n'\
          '136\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.spacegroup == 136
  assert a.shift == 0
  assert a.ifhr == 0
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'Ti' and all(abs(a.atoms[0].pos) < 1e-8)
  assert a.atoms[1].type == 'O' and all(abs(a.atoms[1].pos-[0.30615265, 0.30615265, 0.0]) < 1e-8)
  a.raw = a.raw
  assert a.spacegroup == 136
  assert a.shift == 0
  assert a.ifhr == 0
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'Ti' and all(abs(a.atoms[0].pos) < 1e-8)
  assert a.atoms[1].type == 'O' and all(abs(a.atoms[1].pos-[0.30615265, 0.30615265, 0.0]) < 1e-8)

  a.raw = '1 0 0\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.spacegroup == 'P 42/M N M'
  a.raw = a.raw
  assert a.spacegroup == 136

  a.raw = '1 1 0\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.ifhr == 1
  a.raw = a.raw
  assert a.ifhr == 1

  a.raw = '1 0 1\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.shift == 1
  a.raw = a.raw
  assert a.shift == 1

  a.raw = '1 0 2\n'\
          'P 42/M N M\n'\
          '5 5 5\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert all(abs(a.shift - 5/24.0) < 1e-8)
  a.raw = a.raw
  assert all(abs(a.shift - 5/24.0) < 1e-8)

  print a

if __name__ == "__main__":
  test_crystal()


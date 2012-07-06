def test_general():
  from numpy import all, abs, array
  from quantities import angstrom
  from lada.dftcrystal.basis import Shell
  from lada.error import TypeError
  bhor = 0.5291772083*angstrom # crystal conversion factor
  invbhor = 1e0/bhor/bhor
  invangs = 1e0/angstrom/angstrom

  assert Shell(0).type == 's'  and abs(Shell(0).charge -  2.0) < 1e-8\
         and Shell(0).type_as_int == 0
  assert Shell(1).type == 'sp' and abs(Shell(1).charge -  8.0) < 1e-8\
         and Shell(1).type_as_int == 1
  assert Shell(2).type == 'p'  and abs(Shell(2).charge -  6.0) < 1e-8\
         and Shell(2).type_as_int == 2
  assert Shell(3).type == 'd'  and abs(Shell(3).charge - 10.0) < 1e-8\
         and Shell(3).type_as_int == 3
  assert Shell(4).type == 'f'  and abs(Shell(4).charge - 18.0) < 1e-8\
         and Shell(4).type_as_int == 4

  assert len(Shell(0).append(0.5, 0.1)) == 1
  assert len(Shell(0).append(0.5, 0.1).append(0.6, 0.1)) == 2
  assert abs(Shell(0).append(0.5, 0.1).functions[0][0] - 0.5*invbhor) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][0] - 0.5*invangs) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][0].magnitude - 0.140014258892) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][1] - 0.1) < 1e-8
  assert len(Shell(0).append(0.5*invangs, 0.1).functions[0]) == 2
  assert len(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0]) == 3
  assert abs(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0][1] - 0.1) < 1e-8
  assert abs(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0][2] - 0.2) < 1e-8
  try: Shell(0).append(0.5*invangs, 0.1, 0.2)
  except TypeError: pass
  else: raise Exception()
  try: Shell(1).append(0.5*invangs, 0.1)
  except TypeError: pass
  else: raise Exception()

  a = Shell()
  raw = '0 0 8 2. 1.\n'\
        '225338.0 0.000228\n'\
        '32315.0 0.001929\n'\
        '6883.61 0.011100\n'\
        '1802.14 0.05\n'\
        '543.063 0.17010\n'\
        '187.549 0.369\n'\
        '73.2133 0.4033\n'\
        '30.3718 0.1445\n'
  a.raw = raw
  assert a.type == 's'
  assert abs(a.charge - 2) < 1e-8
  assert len(a) == 8
  assert all(array([len(u) for u in a.functions]) == 2)
  assert all(abs( array([float(u[0]) for u in a.functions])
                  - [225338.0, 32315.0, 6883.61, 1802.14, 543.063, 187.549, 73.2133, 30.3718] ) < 1e-8)
  assert all(abs( array([float(u[1]) for u in a.functions])
                  - [2.28e-4, 1.929e-3, 1.11e-2, 5e-2, 1.701e-1, 3.69e-1, 4.033e-1, 1.445e-1 ] ) < 1e-8)
  a_ = array(a.raw.split(), dtype='float64')
  assert all(abs(array(raw.split(), dtype='float64') - a_) < 1e-8)
  b = eval(repr(a), {'Shell': Shell})
  assert all(abs(array(b.raw.split(), dtype='float64') - a_) < 1e-8)

  raw = '0 1 6 8. 1.\n'\
        '554.042 -0.0059 0.0085\n'\
        '132.525 -0.0683 0.0603\n'\
        '43.6801 -0.1245 0.2124\n'\
        '17.2243 0.2532 0.3902\n'\
        '7.2248 0.6261 0.4097\n'\
        '2.4117 0.282 0.2181\n'
  a.raw = raw
  assert a.type == 'sp'
  assert abs(a.charge - 8) < 1e-8
  assert len(a) == 6
  assert all(array([len(u) for u in a.functions]) == 3)
  assert all(abs( array([float(u[0]) for u in a.functions])
                  - [554.042, 132.525, 43.6801, 17.2243, 7.2248, 2.4117] ) < 1e-8)
  assert all(abs( array([float(u[1]) for u in a.functions])
                  - [-0.0059, -0.0683, -0.1245, 0.2532, 0.6261, 0.282] ) < 1e-8)
  assert all(abs( array([float(u[2]) for u in a.functions])
                  - [0.0085, 0.0603, 0.2124, 0.3902, 0.4097, 0.2181] ) < 1e-8)
  a_ = array(a.raw.split(), dtype='float64')
  assert all(abs(array(raw.split(), dtype='float64') - a_) < 1e-8)
  b = eval(repr(a), {'Shell': Shell})
  assert all(abs(array(b.raw.split(), dtype='float64') - a_) < 1e-8)

if __name__ == "__main__":
  test_general()


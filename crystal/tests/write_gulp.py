def testzb():
  """ Tries and writes a gulp file. """
  from numpy import array, abs, all
  from lada.crystal.binary import zinc_blende
  from lada.crystal.write import gulp

  a = zinc_blende() 
  string = gulp(a).splitlines()
  assert string[0] == 'name'
  assert string[1] == 'Zinc-Blende'
  assert string[2] == 'vectors'
  assert all(abs(array(string[3].split(), dtype='float64') - [0, 0.5, 0.5]) < 1e-8)
  assert all(abs(array(string[4].split(), dtype='float64') - [0.5, 0, 0.5]) < 1e-8)
  assert all(abs(array(string[5].split(), dtype='float64') - [0.5, 0.5, 0]) < 1e-8)
  assert string[6] == 'cartesian'
  assert string[7].split()[:2] == ['A', 'core']
  assert all(abs(array(string[7].split()[2:], dtype='float64')) < 1e-8)
  assert string[8].split()[:2] == ['B', 'core']
  assert all(abs(array(string[8].split()[2:], dtype='float64') - [0.25, 0.25, 0.25]) < 1e-8)

  string2 = gulp(a, symmgroup=216).splitlines()
  assert string2[:6] == string[:6]
  assert string2[6] == 'spacegroup'
  assert string2[7] == '216'
  assert string2[8:] == string[6:]
  
  a.symmgroup = 216
  assert string2 == gulp(a).splitlines()
  
  a[0].asymmetric = True
  a[1].asymmetric = False
  assert string2[:-1] == gulp(a).splitlines()
  del a[0].asymmetric
  assert string2 == gulp(a).splitlines()
  del a[1].asymmetric

  a[1].type = 'A'
  a.symmgroup = 227
  string = gulp(a).splitlines()
  assert string2[:7] == string[:7]
  assert string[7] == '227'
  assert string2[8:-1] == string[8:]

if __name__ == '__main__':
  testzb()

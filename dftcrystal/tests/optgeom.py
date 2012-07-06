def test():
  from collections import namedtuple
  from lada.dftcrystal.optgeom import OptGeom
  Crystal = namedtuple('Crystal', ['optgeom'])

  a = OptGeom()
  assert a.static == True
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.interdun is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg == True
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.interdun is False
  assert a.static is False
  a.maxcycle = None
  assert a.maxcycle is None
  assert a.fulloptg == True
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.interdun is False
  assert a.static is False
  a.cellonly = True
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.interdun is False
  assert a.static is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.interdun is False
  assert a.static is False

  assert a.enabled is False
  assert a.print_input(crystal=Crystal(a)) is None
  a.enabled = True
  assert a.enabled is True
  string = a.print_input(crystal=Crystal(a))
  assert len(string.split('\n')) == 6
  string = string.split('\n')
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = a.print_input(crystal=Crystal(a)).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string

  a.cvolopt = True
  string = a.print_input(crystal=Crystal(a)).split('\n')
  assert len(string) == 7
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = a.print_input(crystal=Crystal(a)).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string
  assert 'CVOLOPT' in string

  a.itatocell = True
  string = a.print_input(crystal=Crystal(a)).split('\n')
  assert len(string) == 6
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = a.print_input(crystal=Crystal(a)).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'ITATOCELL' in string

  a.itatocell = False
  assert a.static
  assert a.print_input() is None
  

if __name__ == "__main__":
  test()


def test():
  from collections import namedtuple
  from lada.dftcrystal.optgeom import OptGeom
  from lada.dftcrystal.input import print_input
  Crystal = namedtuple('Crystal', ['optgeom'])

  a = OptGeom()
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = None
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.cellonly = True
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.intredun is False

  assert a.enabled is False
  assert a.output_map(crystal=Crystal(a)) is None
  a.enabled = True
  assert a.enabled is True
  string = print_input(a.output_map(crystal=Crystal(a)))
  assert len(string.split('\n')) == 6
  string = string.split('\n')
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string

  a.cvolopt = True
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 7
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string
  assert 'CVOLOPT' in string

  a.itatocell = True
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 6
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'ITATOCELL' in string

  a.itatocell = False
  assert a.output_map(crystal=Crystal(a)) is not None
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 5
  assert string[0] == 'OPTGEOM'
  assert string[1] == 'MAXCYCLE'
  assert string[2] == '10'
  assert string[-2] == 'END OPTGEOM'
  assert string[-1] == ''

  a.enabled = False
  assert a.output_map(crystal=Crystal(a)) is None

if __name__ == "__main__":
  test()


def test_fcc():
  from pylada.crystal import bravais
  from pylada.enum import generate_bitstrings
  from fccsets import fccsets

  lattice = bravais.fcc()
  lattice[0].type = ['Si', 'Ge']
  for n in xrange(1, 13):
    result = []
    for x, hft, hermite in generate_bitstrings(lattice, [n]):
      result.append( ''.join(str(i) for i in hermite.flatten()[[0, 3, 4, 6, 7, 8]]) 
                     + ' ' + ''.join(str(i-1) for i in x) )
    assert len(result) == len(fccsets[n])
    assert set(result) == fccsets[n]

def test_ternary():
  from pylada.crystal import bravais
  from pylada.enum import generate_bitstrings
  from ternarysets import ternarysets

  lattice = bravais.fcc()
  lattice[0].type = ['Si', 'Ge', 'C']
  for n in xrange(1, 9):
    result = []
    for x, hft, hermite in generate_bitstrings(lattice, [n]):
      result.append( ''.join(str(i) for i in hermite.flatten()[[0, 3, 4, 6, 7, 8]]) 
                     + ' ' + ''.join(str(i-1) for i in x) )
    assert len(result) == len(ternarysets[n])
    assert set(result) == ternarysets[n]

def test_diamond():
  from pylada.crystal import binary
  from pylada.enum import generate_bitstrings
  from diamondsets import diamondsets

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']
  for n in xrange(2, 8):
    result = []
    for x, hft, hermite in generate_bitstrings(lattice, [n]):
      result.append( ''.join(str(i) for i in hermite.flatten()[[0, 3, 4, 6, 7, 8]]) 
                     + ' ' + ''.join(str(i-1) for i in x) )
    assert len(result) == len(diamondsets[n])
    assert set(result) == diamondsets[n]

def test_zincblende():
  from pylada.crystal import binary
  from pylada.enum import generate_bitstrings
  from zincblendesets import zincblendesets

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Ga']
  for n in xrange(2, 8):
    result = []
    for x, hft, hermite in generate_bitstrings(lattice, [n]):
      result.append( ''.join(str(i) for i in hermite.flatten()[[0, 3, 4, 6, 7, 8]]) 
                     + ' ' + ''.join(str(i-1) for i in x) )
    assert len(result) == len(zincblendesets[n])
    assert set(result) == zincblendesets[n]

if __name__ == '__main__':
  import sys
  sys.path.append(sys.argv[1])
  test_fcc()
  test_ternary()
  test_diamond()
  test_zincblende()

def test():
  from collections import namedtuple
  from lada.dftcrystal.hamiltonian import Dft
  from lada.dftcrystal.parse import parse
  Crystal = namedtuple('Crystal', ['dft'])
  
  a = Dft()
  crystal = Crystal(a)
  assert a.exchange is None
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is False
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert a.print_input(crystal=crystal) is None

  try: a.exchange = 'none'
  except: pass
  else: raise Exception()
  try: a.correlat = 'none'
  except: pass
  else: raise Exception()
  try: a.hybrid = 'none'
  except: pass
  else: raise Exception()
  a.spin = 'none'
  assert a.spin == False

  a.exchange = 'lda'
  assert a.exchange == 'lda'
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is False
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert a.print_input(crystal=crystal) ==  "DFT\nEXCHANGE\nLDA\nEND DFT\n"

  b = parse("DFT\nEXCHANGE\nPBE\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'pbe'
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is False
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc

  a.b3lyp = True
  assert a.exchange == 'becke'
  assert a.correlat == 'lyp'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is False
  assert a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert a.print_input(crystal=Crystal(a)) ==  "DFT\nB3LYP\nEND DFT\n"

  b = parse("DFT\nB3LYP\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.b3lyp
  b = parse("DFT\nPBE0\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.pbe0
  assert a.exchange == 'pbe'
  assert a.correlat == 'pbe'
  assert a.hybrid is None
  assert a.nonlocal.exchange is None
  assert a.nonlocal.correlation is None
  assert a.spin is False
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert a.print_input(crystal=Crystal(a)) ==  "DFT\nPBE0\nEND DFT\n"

  a.b3pw = True
  a.exchange = 'lda'
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is False
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  string = "DFT\nCORRELAT\nPWGGA\nEXCHANGE\nLDA\n"\
           "HYBRID\n" + str(float(20)) + "\n"     \
           "NONLOCAL\n" + str(float(0.9)) + " "   \
           + str(float(0.81)) + "\n" + "END DFT\n"
  assert a.print_input(crystal=Crystal(a)) == string

  b = parse(string, needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is False
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert a.print_input(crystal=Crystal(a)) == string

  string = "DFT\nCORRELAT\nPWGGA\nEXCHANGE\nLDA\n"\
           "HYBRID\n" + str(float(20)) + "\n"     \
           "NONLOCAL\n" + str(float(0.9)) + " "   \
           + str(float(0.81)) + "\nHELLO\nWORLD\n1\n" + "END DFT\n"
  b = parse(string, needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is False
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert 'hello' in a._crysinput
  assert type(a._crysinput['hello']) is bool and a._crysinput['hello'] == True
  assert 'world' in a._crysinput
  assert type(a._crysinput['world']) is int and a._crysinput['world'] == 1
  string = a.print_input(crystal=Crystal(a)).split('\n')
  assert 'HELLO' in string
  assert 'WORLD' in string 
  assert string[string.index('WORLD')+1] == '1'
  a.hello = False
  a.world = 2
  string = a.print_input(crystal=Crystal(a)).split('\n')
  assert 'HELLO' not in string
  assert 'WORLD' in string 
  assert string[string.index('WORLD')+1] == '2'


if __name__ == '__main__': test()

def test():
  from collections import namedtuple
  from pickle import loads, dumps
  from pylada.vasp import Vasp

  Restart = namedtuple('Restart', ['success', 'lmaxmix'])
  a = Vasp()
  o = a._input['lsorbit']
  d = {'LSorbit': o.__class__}
  assert a.lsorbit is None
  assert a.nonscf == False
  assert a._input['lsorbit'].keyword == 'lsorbit'
  assert a._input['nonscf'].keyword is None
  assert o.output_map(vasp=a) is None
  assert eval(repr(o), d).output_map(vasp=a) is None
  assert eval(repr(o), d).value is None
  assert loads(dumps(o)).value is None

  a.lsorbit = True
  assert a.nonscf
  assert a.lsorbit
  try: a._input['lsorbit'].output_map(vasp=a)
  except ValueError: pass
  else: raise Exception()
  a.restart = Restart(False, 7)
  try: a._input['lsorbit'].output_map(vasp=a)
  except ValueError: pass
  else: raise Exception()
  a.restart = Restart(True, 7)
  assert 'lsorbit' in o.output_map(vasp=a)
  assert o.output_map(vasp=a)['lsorbit'] == '.TRUE.'
  assert a.lmaxmix == 7
  a.lmaxmix = 5
  a.restart = Restart(True, 6)
  assert 'lsorbit' in o.output_map(vasp=a)
  assert o.output_map(vasp=a)['lsorbit'] == '.TRUE.'
  assert a.lmaxmix == 6
  assert loads(dumps(o)).value is True
  assert eval(repr(o), d).value is True
  

if __name__ == "__main__": test()
  

def test():
  from pickle import loads, dumps
  from lada.vasp import Vasp
  from lada.crystal import Structure
  a = Vasp()


  o = a._input['ediff']
  d = {'Ediff': o.__class__}
  assert a.ediff is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'ediff'
  assert loads(dumps(o)).output_map() is None

  u = 0.25
  x, y = u, 0.25-u
  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
                       .add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "Mg") \
                       .add_atom(5.000000e-01, 2.500000e-01, 2.500000e-01, "Mg") \
                       .add_atom(2.500000e-01, 5.000000e-01, 2.500000e-01, "Mg") \
                       .add_atom(2.500000e-01, 2.500000e-01, 5.000000e-01, "Mg") \
                       .add_atom(8.750000e-01, 8.750000e-01, 8.750000e-01, "Al") \
                       .add_atom(1.250000e-01, 1.250000e-01, 1.250000e-01, "Al") \
                       .add_atom(     x,     x,     x, "O") \
                       .add_atom(     x,     y,     y, "O") \
                       .add_atom(     y,     x,     y, "O") \
                       .add_atom(     y,     y,     x, "O") \
                       .add_atom(    -x,    -x,    -x, "O") \
                       .add_atom(    -x,    -y,    -y, "O") \
                       .add_atom(    -y,    -x,    -y, "O") \
                       .add_atom(    -y,    -y,    -x, "O") 
  
  N = float(len(structure))
  a.ediff = 2e-4
  assert N < 100
  assert abs(a.ediff - 2e-4) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediff']) - a.ediff * N) < a.ediff * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediff']) - a.ediff * N) < a.ediff * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediff']) - a.ediff * N) < a.ediff * 1e-2
  a.ediff = -2e-4
  assert N < 100
  assert abs(a.ediff + 2e-4) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediff']) + a.ediff) < -a.ediff * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediff']) + a.ediff) < -a.ediff * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediff']) + a.ediff) < -a.ediff * 1e-2

  a.ediffg = -2e-4
  o = a._input['ediffg']
  d = {'Ediffg': o.__class__}
  assert N < 100
  assert abs(a.ediffg + 2e-4) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2
 
if __name__ == '__main__': test()

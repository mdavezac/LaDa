def test_ediff():
  from pickle import loads, dumps
  from lada.vasp import Vasp
  from lada.crystal import Structure
  a = Vasp()

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

  o = a._input['ediff']
  d = {'Ediff': o.__class__}
  assert a.ediff is None
  assert a.ediff_per_atom is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'ediff'
  assert loads(dumps(o)).output_map() is None
  a.ediff_per_atom = 1e-5
  a.ediff = 2e-4
  assert abs(a.ediff - 2e-4) < 1e-8
  assert a.ediff_per_atom is None
  assert abs(float(o.output_map(structure=structure)['ediff']) - a.ediff) < a.ediff * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediff']) - a.ediff) < a.ediff * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediff']) - a.ediff) < a.ediff * 1e-2
  a.ediff = -1
  assert abs(a.ediff) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediff'])) < 1e-8

  a = Vasp()
  o = a._input['ediff_per_atom']
  d = {'EdiffPerAtom': o.__class__}
  assert a.ediff_per_atom is None
  assert a.ediff is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'ediff'
  assert loads(dumps(o)).output_map() is None
  a.ediff = 1e-5
  a.ediff_per_atom = 2e-4
  assert abs(a.ediff_per_atom - 2e-4) < 1e-8
  assert a.ediff is None
  assert abs(float(o.output_map(structure=structure)['ediff']) - 2e-4 * N) < 2e-4 * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediff']) - 2e-4 * N) < 2e-4 * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediff']) - 2e-4 * N) < 2e-4 * 1e-2

  a.ediff = 1e-4
  a.ediff_per_atom = None
  assert abs(a.ediff-1e-4) < 1e-8
  assert a.ediff_per_atom is None
  a.ediff_per_atom = 1e-4
  a.ediff = None
  assert abs(a.ediff_per_atom-1e-4) < 1e-8
  assert a.ediff is None


def test_ediffg():
  from pickle import loads, dumps
  from lada.vasp import Vasp
  from lada.crystal import Structure
  a = Vasp()

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

  o = a._input['ediffg']
  d = {'Ediffg': o.__class__}
  assert a.ediffg is None
  assert a.ediffg_per_atom is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'ediffg'
  assert loads(dumps(o)).output_map() is None
  a.ediffg_per_atom = 1e-5
  a.ediffg = 2e-4
  assert abs(a.ediffg - 2e-4) < 1e-8
  assert a.ediffg_per_atom is None
  assert abs(float(o.output_map(structure=structure)['ediffg']) - a.ediffg) < a.ediffg * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediffg']) - a.ediffg) < a.ediffg * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediffg']) - a.ediffg) < a.ediffg * 1e-2
  a.ediffg_per_atom = 1e-5
  a.ediffg = -2e-4
  assert abs(a.ediffg + 2e-4) < 1e-8
  assert a.ediffg_per_atom is None
  assert abs(float(o.output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediffg']) - a.ediffg) < -a.ediffg * 1e-2


  a.ediffg = 2e-4
  a.ediffg_per_atom = 2e-4
  o = a._input['ediffg_per_atom']
  d = {'EdiffgPerAtom': o.__class__}
  assert a.ediffg is None
  assert abs(a.ediffg_per_atom - 2e-4) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediffg']) - 2e-4*N) < 2e-4 * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediffg']) - 2e-4*N) < 2e-4 * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediffg']) - 2e-4*N) < 2e-4 * 1e-2
  a.ediffg = 2e-4
  a.ediffg_per_atom = -2e-4
  assert a.ediffg is None
  assert abs(a.ediffg_per_atom + 2e-4) < 1e-8
  assert abs(float(o.output_map(structure=structure)['ediffg']) + 2e-4) < 2e-4 * 1e-2
  assert abs(float(eval(repr(o), d).output_map(structure=structure)['ediffg']) + 2e-4) < 2e-4 * 1e-2
  assert abs(float(loads(dumps(o)).output_map(structure=structure)['ediffg']) + 2e-4) < 2e-4 * 1e-2
  a.ediffg_per_atom = None
  assert a.ediffg is None
  assert a.ediffg_per_atom is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'ediffg'
  assert loads(dumps(o)).output_map() is None

  a.ediffg = 1e-4
  a.ediffg_per_atom = None
  assert abs(a.ediffg-1e-4) < 1e-8
  assert a.ediffg_per_atom is None
  a.ediffg_per_atom = 1e-4
  a.ediffg = None
  assert abs(a.ediffg_per_atom-1e-4) < 1e-8
  assert a.ediffg is None

if __name__ == '__main__':
  test_ediff()
  test_ediffg()

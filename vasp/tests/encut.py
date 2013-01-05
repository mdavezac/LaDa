def test(path):
  from pickle import loads, dumps
  from quantities import eV, hartree
  import quantities 
  import numpy
  from pylada.vasp import Vasp
  from pylada.crystal import Structure

  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")
  a = Vasp()
  a.add_specie = "Si", "{0}/pseudos/Si".format(path)


  o = a._input['encut']
  d = {'Encut': o.__class__}
  d.update(quantities.__dict__)
  d.update(numpy.__dict__)
  assert a.ediff is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'encut'
  assert loads(dumps(o)).output_map() is None

  a.encut = 1e0
  assert abs(a.encut - 1e0) < 1e-8
  assert abs(float(o.output_map(structure=structure, vasp=a)['encut'])-245.345) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encut'])-245.345) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encut'])-245.345) < 1e-8
  assert abs(eval(repr(o), d).value - 1.0) < 1e-8
  assert abs(loads(dumps(o)).value - 1.0) < 1e-8
  a.encut = 0.8
  assert abs(a.encut - 0.8) < 1e-8
  assert abs(float(o.output_map(structure=structure, vasp=a)['encut'])-245.345*0.8) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encut'])-245.345*0.8) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encut'])-245.345*0.8) < 1e-8
  assert abs(eval(repr(o), d).value - 0.8) < 1e-8
  assert abs(loads(dumps(o)).value - 0.8) < 1e-8
  a.encut = 200
  assert abs(a.encut - 200) < 1e-8
  assert abs(float(o.output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(eval(repr(o), d).value - 200) < 1e-8
  assert abs(loads(dumps(o)).value - 200) < 1e-8
  a.encut = 200 * eV
  assert abs(a.encut - 200 * eV) < 1e-8
  assert a.encut.units == eV
  assert abs(float(o.output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(eval(repr(o), d).value - 200 * eV) < 1e-8
  assert abs(loads(dumps(o)).value - 200 * eV) < 1e-8
  assert eval(repr(o), d).value.units == eV
  assert loads(dumps(o)).value.units == eV
  a.encut = (200 * eV).rescale(hartree)
  assert a.encut.units == hartree
  assert abs(a.encut - 200 * eV) < 1e-8
  assert abs(float(o.output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encut'])-200) < 1e-8
  assert abs(eval(repr(o), d).value - 200 * eV) < 1e-8
  assert abs(loads(dumps(o)).value - 200 * eV) < 1e-8
  assert eval(repr(o), d).value.units == hartree
  assert loads(dumps(o)).value.units == hartree

  o = a._input['encutgw']
  d = {'EncutGW': o.__class__}
  d.update(quantities.__dict__)
  d.update(numpy.__dict__)
  assert a.ediff is None
  assert o.output_map() is None
  assert eval(repr(o), d).output_map() is None
  assert eval(repr(o), d).keyword == 'encutgw'
  assert loads(dumps(o)).output_map() is None
  a.encutgw = (200 * eV).rescale(hartree)
  assert a.encutgw.units == hartree
  assert abs(a.encutgw - 200 * eV) < 1e-8
  assert abs(float(o.output_map(structure=structure, vasp=a)['encutgw'])-200) < 1e-8
  assert abs(float(eval(repr(o), d).output_map(structure=structure, vasp=a)['encutgw'])-200) < 1e-8
  assert abs(float(loads(dumps(o)).output_map(structure=structure, vasp=a)['encutgw'])-200) < 1e-8
  assert abs(eval(repr(o), d).value - 200 * eV) < 1e-8
  assert abs(loads(dumps(o)).value - 200 * eV) < 1e-8
  assert eval(repr(o), d).value.units == hartree
  assert loads(dumps(o)).value.units == hartree
 
if __name__ == "__main__":
  from sys import argv
  test(argv[1])
  

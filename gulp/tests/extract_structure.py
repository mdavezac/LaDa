def test_noopti_anatase(path):
  """ Test structure extraction. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom
  from lada.gulp import Extract

  extract = Extract(join(path, 'noopti.gout'))
  assert not extract.optimize
  assert not extract.conp
  assert extract.qeq
  cell = array([[-1.897114,  1.897114,  4.884209],
                [ 1.897114, -1.897114,  4.884209],
                [ 1.897114,  1.897114, -4.884209]]) * angstrom
  assert all(abs(extract.cell - cell) < 1e-4)
  assert extract.cell is extract.input_cell
  assert all(abs(extract.structure.cell - cell.magnitude) < 1e-4)
  assert len(extract.structure) == 2
  assert extract.structure[0].type == 'Ti'
  assert all(abs(extract.structure[0].pos) < 1e-4)
  assert extract.structure[1].type == 'O'
  assert all(abs(extract.structure[1].pos - [0, 0, 0.20474]) < 1e-4)
  assert all(abs(extract.structure[1].pos - [0, 0, 0.20474]) < 1e-4)
  assert extract.structure is extract.input_structure
  assert getattr(extract.structure[0], 'label', 0) == 1
  assert getattr(extract.structure[1], 'label', 0) == 2
  assert abs(getattr(extract.structure[0], 'occupancy', 0) - 1e0) < 1e-8
  assert abs(getattr(extract.structure[1], 'occupancy', 0) - 1e0) < 1e-8


def test_opti_anatase(path):
  """ Test structure optimization extraction. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom, e, eV
  from lada.gulp import Extract

  extract = Extract(join(path, 'opti.gout'))
  assert extract.optimize
  assert extract.conp
  assert extract.qeq
  cell = array([[-1.897114,  1.897114,  4.884209],
                [ 1.897114, -1.897114,  4.884209],
                [ 1.897114,  1.897114, -4.884209]]) * angstrom
  assert all(abs(extract.input_cell - cell) < 1e-4)
  assert extract.cell is not extract.input_cell
  assert all(abs(extract.input_structure.cell - cell.magnitude) < 1e-4)
  assert len(extract.input_structure) == 2
  assert extract.input_structure[0].type == 'Ti'
  assert all(abs(extract.input_structure[0].pos) < 1e-4)
  assert extract.input_structure[1].type == 'O'
  assert all(abs(extract.input_structure[1].pos - [0, 0, 0.20474]) < 1e-4)
  assert all(abs(extract.input_structure[1].pos - [0, 0, 0.20474]) < 1e-4)
  assert extract.structure is not extract.input_structure

  cell = array([[-1.924847,  1.924847,  4.53005 ],
                [ 1.924847, -1.924847,  4.53005 ],
                [ 1.924847,  1.924847, -4.53005 ]]) * angstrom
  assert all(abs(extract.cell - cell) < 1e-4)
  assert len(extract.structure) == 2
  assert extract.structure[0].type == 'Ti'
  assert all(abs(extract.structure[0].pos) < 1e-4)
  assert extract.structure[1].type == 'O'
  assert all(abs(extract.structure[1].pos - [0, 0, 0.216931]) < 1e-4)
  assert getattr(extract.structure[0], 'label', 0) == 1
  assert getattr(extract.structure[1], 'label', 0) == 2
  assert abs(getattr(extract.structure[0], 'charge', 0) - 1.131064*e) < 1e-5
  assert abs(getattr(extract.structure[1], 'charge', 0) + 0.565532*e) < 1e-5
  assert abs(extract.structure.energy + 16.38287598 * eV) < 1e-5

if __name__ == '__main__':
  from sys import argv
  from os.path import dirname, join
  test_noopti_anatase(join(dirname(argv[0]), 'data'))
  test_opti_anatase(join(dirname(argv[0]), 'data'))

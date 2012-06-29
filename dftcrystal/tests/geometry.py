def test_crystal():
  from numpy import array
  from lada.dftcrystal.geometry import Crystal
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O')
  assert structure.print_input() is not None
  assert len(structure.print_input().split('\n')) == 9
  assert all(abs(array(structure.print_input().split()[1:-2], dtype='float64') \
             - [0, 0, 0, 136, 4.63909875, 2.97938395, 2, 22, 0, 0, 0, 8, 0.306153, 0.306153, 0]) < 1e-8)
  assert structure.print_input().split()[0] == 'CRYSTAL'
  assert structure.print_input().split()[-2] == 'END'
  assert structure.print_input().split()[-1] == 'CRYSTAL'

if __name__ == '__main__':
  test_crystal()

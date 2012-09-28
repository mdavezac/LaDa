def test_crystal():
  from numpy import array
  from lada.dftcrystal.crystal import Crystal
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

def test_supercell():
  from numpy import array
  from lada.dftcrystal import Supercell, Crystal
  from lada.error import ValueError
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O') \
                  .append(Supercell([[2,0,0],[0,1,0],[0,0,1]]))
  assert len(structure.eval()) == 2 * 6
  structure[-1].matrix[2,2] = 2
  assert len(structure.eval()) == 2 * 2 * 6

  # test errors
  try: structure[-1].matrix = array( [[0]*3]*3 )
  except ValueError: pass
  else: raise Exception()
  structure[-1].matrix = array([[1,0],[0,1]])
  try: structure.print_input()
  except ValueError: pass
  else: raise Exception()



if __name__ == '__main__':
  test_crystal()
  test_supercell()

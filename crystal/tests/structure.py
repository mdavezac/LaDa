""" Checks structure methods and attributes. """
def test_init(Class):
  """ Test structure initialization. """
  from numpy import all, abs, array, identity

  a = Class()
  assert all(abs(a.cell - identity(3)) < 1e-8) and abs(a.scale - 1e0) < 1e0\
         and len(a.__dict__) == 0

  a = Class(identity(3)*2.5, scale=5.45)
  assert all(abs(a.cell - identity(3)*2.5) < 1e-8) and abs(a.scale - 5.45) < 1e0\
         and len(a.__dict__) == 0

  a = Class(identity(3)*2.5, scale=5.45, m=True)
  assert all(abs(a.cell - identity(3)*2.5) < 1e-8) and abs(a.scale - 5.45) < 1e0\
         and len(a.__dict__) == 1 and getattr(a, 'm', False)
  assert str(a) == "structure = Structure( 2.5, 0, 0,\\\n"\
                   "                       0, 2.5, 0,\\\n"\
                   "                       0, 0, 2.5,\\\n"\
                   "                       scale=5.45, m=True )"
  a.add_atom(0,0,0, "Au")\
            (0.25, 0.5, 0.25, "Au", "Pd", m=5)
  assert str(a) == "structure = Structure( 2.5, 0, 0,\\\n"\
                   "                       0, 2.5, 0,\\\n"\
                   "                       0, 0, 2.5,\\\n"\
                   "                       scale=5.45, m=True )\n"\
                   "structure.add_atom(0, 0, 0, 'Au')\\\n"\
                   "                  (0.25, 0.5, 0.25, ['Au', 'Pd'], m=5)"
  a.cell[0,0] = 1e0
  a.cell[1,:] = 1e0
  assert str(a) == "structure = Structure( 1, 1, 0,\\\n"\
                   "                       0, 1, 0,\\\n"\
                   "                       0, 1, 2.5,\\\n"\
                   "                       scale=5.45, m=True )\n"\
                   "structure.add_atom(0, 0, 0, 'Au')\\\n"\
                   "                  (0.25, 0.5, 0.25, ['Au', 'Pd'], m=5)"


if __name__ == "__main__":
  from lada.crystal.cppwrappers import Structure
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  # tries to run test with normal class.
  test_init(Structure) 

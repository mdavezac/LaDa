""" Checks split-config routine. """
def check_bcc():
  from numpy import abs
  from pylada.crystal.cppwrappers import Structure, splitconfigs, supercell
  structure = Structure([[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]])\
                       .add_atom(0,0,0,"Mo")
  configs = splitconfigs(structure, structure[0], 12)
  assert len(configs) == 1
  assert abs(configs[0][1] - 1e0) < 1e-8
  for u in configs[0][0]: assert u[0] is structure[0]
  assert all(abs(configs[0][0][0][1] - [ 0.,  0.,  0.]) < 1e-8)
  assert all(abs(configs[0][0][1][1] - [  8.66025404e-01,  -1.11022302e-16,   0.00000000e+00]) < 1e-8)
  assert all(abs(configs[0][0][2][1] - [ 0.28867513,  0.81649658,  0.        ]) < 1e-8)
  assert all(abs(configs[0][0][3][1] - [ 0.28867513, -0.40824829,  0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][4][1] - [ 0.28867513, -0.40824829, -0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][5][1] - [-0.28867513,  0.40824829,  0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][6][1] - [-0.28867513,  0.40824829, -0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][7][1] - [-0.28867513, -0.81649658,  0.        ]) < 1e-8)
  assert all(abs(configs[0][0][8][1] - [ -8.66025404e-01,   1.11022302e-16,   0.00000000e+00]) < 1e-8)
  assert all(abs(configs[0][0][9][1] - [ 0.57735027,  0.40824829,  0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][10][1] - [ 0.57735027,  0.40824829, -0.70710678]) < 1e-8)
  assert all(abs(configs[0][0][11][1] - [ 0.57735027, -0.81649658,  0.        ]) < 1e-8)

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  check_bcc()

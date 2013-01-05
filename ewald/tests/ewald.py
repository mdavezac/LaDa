def test():
  from numpy import all, abs, sqrt
  from pylada.crystal import Structure
  from pylada.ewald import ewald
  from quantities import angstrom, a0, Ry
  structure = Structure( [ [1,0,0],
                           [0,1,0],
                           [0,0,1] ], scale=50 )                               \
                       .add_atom(0, 0, 0, 'A', charge=1e0)                     \
                       .add_atom( float(a0.rescale(angstrom)/50.0), 0, 0, 'A',  
                                  charge=-1e0 )
  result = ewald(structure, cutoff = 80)
  assert abs(result.energy + 2e0*Ry) < 1e-3
  assert all(abs(abs(result[0].force) - [2e0, 0, 0]*Ry/a0) < 1e-3)
  assert all(abs(abs(result[1].force) - [2e0, 0, 0]*Ry/a0) < 1e-3)

  a = float(a0.rescale(angstrom)/50.0) / sqrt(2.)
  structure = Structure( [ [1,0,0],
                           [0,1,0],
                           [0,0,1] ], scale=50 )                               \
                       .add_atom(0, 0, 0, 'A', charge=1e0)                     \
                       .add_atom(0, a, a, 'A', charge=-1e0 )
  result = ewald(structure, cutoff = 80)
  assert abs(result.energy + 2e0*Ry) < 1e-3
  assert all(abs(abs(result[0].force) - [0, 2./sqrt(2), 2./sqrt(2)]*Ry/a0) < 1e-3)
  assert all(abs(abs(result[1].force) - [0, 2./sqrt(2), 2./sqrt(2)]*Ry/a0) < 1e-3)

if __name__ == '__main__':
  test()

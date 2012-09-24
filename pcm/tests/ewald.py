def test():
  import lada
  from numpy import all, abs
  from quantities import angstrom, a0, Ry
  from lada.crystal import Structure
  from lada.pcm import ewald
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

if __name__ == '__main__':
  test()

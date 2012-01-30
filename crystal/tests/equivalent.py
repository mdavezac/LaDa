""" Checks crystal method `lada.crystal.equivalent`. """
def scale(A, B):
  """ Check changes to scale. """
  from lada.crystal.cppwrappers import equivalent
  B = B.copy()
  assert equivalent(A, B)
  B.scale = 3.0;
  assert not equivalent(A, B)
  assert equivalent(A, B, scale=False)
  B.cell *= 0.5;
  for atom in B: atom.pos *= 0.5
  assert not equivalent(A, B)
  assert equivalent(A, B, scale=False)

def motif(A, B):
  """ Check changes in motif. """
  from numpy import dot
  from lada.crystal.cppwrappers import cell_invariants
  for op in cell_invariants(A):
    B = B.copy()
    for atom in B:
      atom.pos = dot(op[:3] * atom.pos)
      scale(A, B)

def basis(A, B):
  """ Adds rotation and translation of cartesian basis. """
  from numpy import dot, pi
  from lada.crystal.cppwrappers import cell_invariants, transform
  from lada.math import Translation, Rotation
  from random import random
  motif(A, transform(B, Rotation(0.5 * pi, [1,0,0])))
  motif(A, transform(B, Rotation(-pi, [1,0,0])))
  motif(A, transform(B, Rotation(-0.13*pi, [1,0,0])))
  motif(A, transform(B, Translation([0.25, 0.25, 0.25])))
  motif(A, transform(B, Rotation(random()*2*pi, [1, 0, 0]) \
                        + Translation([random()-0.5,random()-0.5,random()-0.5])))
def decoration(A, B, lattice):
  """ Adds changes to the motif. """
  from lada.crystal.cppwrappers import SmithTransform
  from lada.math import is_integer
  smith = SmithTransform(lattice, A)

  # create map of atoms.
  indices = [-1] * len(A)
  for atom, i in enumerate(A):
    indices[smith.index(atom.pos-lattice[atom.site].pos, atom.site)] = i

  # transform A according to all possible atom-atom translations.
  for atom in A:
    if atom.site != A[0].site: continue # only primitive translations.
    trans = atom.pos - A[0].pos
    B = A.copy()
    for index in indices:
      vec = A[index].pos + trans - lattice[A[index].site].pos
      assert is_integer(dot(inv(lattice.cell), vec)) # vector should be on lattice
      B[ indices[smith.index(vec, A[index].site)] ] = A[index]
    basis(A, B)

def test():
  from lada.crystal.cppwrappers import Structure

  zb = Structure( 0,0.5,0.5,
                  0.5,0,0.5,
                  0.5,0.5,0 )\
                .add_atom(0,0,0, "Si")\
                .add_atom(0,0,0, "Si", "Ge")
  basis(zb, zb)

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

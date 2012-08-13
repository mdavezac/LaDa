""" Checks that space group is correct. """
def test_fcc():
  """ Test fcc space-group and equivalents """
  from numpy import all, abs, dot
  from lada.crystal.cppwrappers import space_group, Structure, equivalent, transform

  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]], m=True).add_atom(0,0,0,"Si", m=True)
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)
    assert all(abs(op[3, :]) < 1e-8) 

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    assert getattr(other, 'm', False) 
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type
      assert getattr(atom, 'm', False) 

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)

def test_b5(u):
  """ Test b5 space-group and equivalents """
  from numpy import all, abs, dot
  from lada.crystal.cppwrappers import space_group, Structure, equivalent, transform

  x, y = u, 0.25-u
  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
                       .add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "A") \
                       .add_atom(5.000000e-01, 2.500000e-01, 2.500000e-01, "A") \
                       .add_atom(2.500000e-01, 5.000000e-01, 2.500000e-01, "A") \
                       .add_atom(2.500000e-01, 2.500000e-01, 5.000000e-01, "A") \
                       .add_atom(8.750000e-01, 8.750000e-01, 8.750000e-01, "B") \
                       .add_atom(1.250000e-01, 1.250000e-01, 1.250000e-01, "B") \
                       .add_atom(     x,     x,     x, "X") \
                       .add_atom(     x,     y,     y, "X") \
                       .add_atom(     y,     x,     y, "X") \
                       .add_atom(     y,     y,     x, "X") \
                       .add_atom(    -x,    -x,    -x, "X") \
                       .add_atom(    -x,    -y,    -y, "X") \
                       .add_atom(    -y,    -x,    -y, "X") \
                       .add_atom(    -y,    -y,    -x, "X") 
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)

  structure[0], structure[-1] = structure[-1], structure[0]
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)
    
def test_zb():
  from numpy import all, abs, dot
  from lada.crystal import space_group, transform, binary
  from lada.crystal.cppwrappers import equivalent

  structure = binary.zinc_blende()
  ops = space_group(structure)
  assert len(ops) == 24
  for op in ops:
    assert op.shape == (4, 3)

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
    for a, atom in zip(structure, other):
      assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
      assert a.type == atom.type

    assert equivalent(structure, other, cartesian=False)
    assert equivalent(other, structure, cartesian=False)
     
  for atom in structure: atom.type = ['A', 'B']
  ops = space_group(structure)
  assert len(ops) == 48
  for op in ops:
    assert op.shape == (4, 3)

    other = transform(structure, op)
    assert all(abs(dot(op[:3], structure.cell)-other.cell) < 1e-8)
 #  for a, atom in zip(structure, other):
 #    assert all(abs(dot(op[:3], a.pos) + op[3] - atom.pos) < 1e-8)
 #    assert a.type == atom.type

 #  assert equivalent(structure, other, cartesian=False)
 #  assert equivalent(other, structure, cartesian=False)


if __name__ == "__main__":
  test_fcc()
  test_b5(0.25)
  test_b5(0.36)
  test_zb()

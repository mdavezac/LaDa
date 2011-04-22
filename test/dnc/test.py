from numpy import diag, array
from numpy.linalg import norm
from matplotlib.pyplot import scatter, show
from lada.crystal.binary import zinc_blende
from lada.crystal import dnc_iterator

structure = zinc_blende().to_structure(diag([20,10,10]))
result = [False for u in structure.atoms]
for box in dnc_iterator(structure, 15, 0.3):

  all = [u for u in box]
  small_box = [structure.atoms[i].pos + trans for i, trans, f in all if f == True]
  large_box = array([structure.atoms[i].pos + trans for i, trans, f in all])
  for index, trans, in_small_box in all:
    assert result[index] == False
    result[index] = False

  for pos in small_box:
    sizes = sorted([ norm(u) for u in (large_box - pos) ])
    assert sizes.pop(0) < 1e-12
    for i in range(1, 4): 
      assert abs(sizes[i] - sizes[0]) < 1e-12
    assert abs(sizes[4] - sizes[0]) > 1e-12
    
    

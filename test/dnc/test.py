from numpy import diag, array, dot
from numpy.linalg import norm, inv
from matplotlib.pyplot import scatter, show
from numpy import array
from operator import itemgetter
from lada.crystal.binary import zinc_blende
from lada.crystal import dnc_iterator

cell = array([[10.0, 0.5, 0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0]], dtype="float64")
structure = zinc_blende().to_structure(cell) #diag([10,1,1]))
# print structure
# print dnc_iterator(structure, 45, 1.5).mesh
result = [False for u in structure.atoms]
for box in dnc_iterator(structure, 15, 0.9):

  all = [u for u in box]
  small_box = [(j, structure.atoms[i].pos + trans) for j, (i, trans, f) in enumerate(all) if f == True]
  large_box = array([structure.atoms[i].pos + trans for i, trans, f in all])
  for index, trans, in_small_box in all:
    assert result[index] == False
    result[index] = False

  for index, pos in small_box:
    sizes = sorted([ (i, norm(u)) for i, u in enumerate(large_box - pos) ], key=itemgetter(1))
#   print structure.atoms[all[index][0]].pos, all[index][0]
#   for i, bondlength in sizes:
#     print " ", bondlength, structure.atoms[all[i][0]].pos, structure.atoms[all[i][0]].pos + all[i][1]
    for i in range(2, 5): 
      assert abs(sizes[i][1] - sizes[1][1]) < 1e-12
    assert abs(sizes[5][1] - sizes[1][1]) > 1e-12
    
    

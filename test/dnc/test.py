from numpy import diag, array
from matplotlib.pyplot import scatter, show
from lada.crystal.binary import zinc_blende
from lada.crystal import dnc_iterator

structure = zinc_blende().to_structure(diag([20,10,10]))
colors = "rgbcmyk"
scatter( [a.pos[0] for a in structure.atoms], [a.pos[1] for a in structure.atoms], s=200, linewidth=0 )
ibox = dnc_iterator(structure, 15, 0.3)
box = ibox.next()
box = ibox.next()
box = ibox.next()
box = ibox.next()
result = [], []
for index, translation, in_small_box in box:
  if in_small_box: 
    result[0].append(structure.atoms[index].pos+translation)
  else:
    print translation
    result[1].append(structure.atoms[index].pos+translation)

all = []
for box in  dnc_iterator(structure, 15, 0.5):
  for index, translation, in_small_box in box:
    if not in_small_box: all.append(structure.atoms[index].pos + translation)
# scatter( [a[0] for a in all], [a[1] for a in all],  c='g', s=120, linewidth=0)
scatter( [a[0] for a in result[1]], [a[1] for a in result[1]],  c='g', s=120, linewidth=0)
scatter( [a[0] for a in result[0]], [a[1] for a in result[0]],  c='r', s=50, linewidth=0)
show()

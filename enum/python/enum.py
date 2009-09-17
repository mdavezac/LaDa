#! /usr/bin/python

def add_directories():
  import sys
  import os

  i = 0
# for i, p in enumerate(sys.path):
#   if p == '/home/davezac/usr/lib64/python2.5/site-packages':
#     break;
# sys.path.pop(i)

  sys.path.extend\
  (\
    [\
     '/home/davezac/usr/debug/lib64/python2.5/site-packages'\
    ]\
  )

def find_smith_set( _n, _hermite ):
  from math import pow
  file = open("struct_enum.out", "r") 
  result = []
  for line in file:
    if line.split()[0] == "#tot": break
  for line in file:
    splitted = line.split()
    if int(splitted[2]) != _n: continue
    if int( splitted[7]) != _hermite[0,0]: continue
    if int( splitted[8]) != _hermite[1,0]: continue
    if int( splitted[9]) != _hermite[1,1]: continue
    if int(splitted[10]) != _hermite[2,0]: continue
    if int(splitted[11]) != _hermite[2,1]: continue
    if int(splitted[12]) != _hermite[2,2]: continue
   
    number = 0
    base = int(pow(2, _hermite[0,0]*_hermite[1,1]*_hermite[2,2]-1))
    for character in splitted[-1]:
      number += int(character) * base
      base /= 2
    result.append(number)
  return result


def create_lattice():

  from lada import crystal, atat
  from math import sqrt

  lattice = crystal.Lattice()

  lattice.cell = atat.rMatrix3d( [ [           0.5,          0.5,         0 ], \
                                   [  -sqrt(3)/2.0,  sqrt(3)/2.0,         0 ], \
                                   [             0,            0,         1 ] ] )
  lattice.scale = 4.42

  lattice.sites.append( crystal.Site( (0, 0,0) ) )
  lattice.sites[0].type = crystal.StringVector( [ "K", "Rb" ] );
  lattice.sites.append( crystal.Site( (0, -2.0/sqrt(12), 0.5) ) )
  lattice.sites[1].type = crystal.StringVector( [ "K", "Rb" ] );
  lattice.find_space_group()

  return lattice

def all_symmetrics( _x, _trans, _label, _rotations, _fl ):

  def all_smith( __x):
    for labelperm in _label: yield labelperm(__x, _fl), 'l'
    for translation in _trans:
      u = translation(__x, _fl)
      yield u, 't'
      for labelperm in _label: yield labelperm(u, _fl), 'l'

  for u in all_smith(_x): yield u
  for transform in _rotations:
    t = transform(_x, _fl)
    yield t, 'r'
    for u in all_smith(t): yield u


# def main():
from lada import enumeration, atat
from math import pow
import time

lattice = create_lattice()
lattice.set_as_crystal_lattice()
species = ("K", "Rb")

nconf = 1
t0 = time.time()
for n in range(2, 5):
  supercells = enumeration.find_all_cells(lattice, n)
  smiths = enumeration.create_smith_groups(lattice, supercells)
  nflavors = enumeration.count_flavors( lattice )
  nsites = len(lattice.sites)
  transforms = enumeration.create_transforms(lattice)
  for smith in smiths:
    card = smith.smith[0]*smith.smith[1]*smith.smith[2]*nsites
    label_exchange=enumeration.LabelExchange( card, nflavors )
    flavorbase = enumeration.create_flavorbase(card, nflavors)
    translations = enumeration.Translation(smith.smith, nsites)
    database = enumeration.Database(card, nflavors)
    maxterm = 0
    for x in xrange(1, int(pow(nflavors, card))-1):
      if not database[x]: continue
      maxterm = x
      for labelperm in label_exchange:
        t = labelperm(x, flavorbase)
        if t > x: database[t] = False
      for translation in translations:
        t = translation(x, flavorbase)
        if t > x: database[t] = False
        elif t == x: 
          database[t] = False
          continue
        for labelperm in label_exchange:
          u = labelperm(t, flavorbase)
          if u > x: database[u] = False

    # checks supercell dependent transforms.
    for nsupercell, supercell in enumerate(smith.supercells):
      mine = []
      # creates list of transformation which leave the supercell invariant.
      cell = lattice.cell * supercell.hermite
      specialized = []
      for transform in transforms:
        if transform.trans != atat.rVector3d(0,0,0): print transform.trans
        if not transform.invariant(cell): continue
        transform.init(supercell.transform, smith.smith)
        if not transform.is_trivial:
          specialized.append( transform )

      specialized_database = enumeration.Database(database)
      for x in xrange(1, maxterm+1):
        if not database[x]: continue
        maxterm = x
        
        for transform in specialized:
          t = transform(x, flavorbase)
          if t == x: continue
          specialized_database[t] = False

          for labelperm in label_exchange:
            u = labelperm(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            
          for translation in translations:
            u = translation(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            for labelperm in label_exchange:
              v = labelperm(u, flavorbase)
              if v == x: continue
              specialized_database[v] = False
        if specialized_database[x]: # mine.append(x)
          print "%5i %2i " % (nconf, n),
          for i in range(0,3):
            for j in range(0,i+1):
              print "%2i " % (supercell.hermite[i,j]),
          print "%10i %20s " % (int(x), enumeration.as_bitstring(x, flavorbase))
        
#     gus = set(find_smith_set(n, supercell.hermite));
#     mine = set(mine)
#     print "smith: ", smith.smith[0], smith.smith[1], smith.smith[1], "supercell: ", supercell, \
#           maxterm
#     for x in gus-mine:
#       is_not_there = True
#       for u, type in all_symmetrics(x, translations, label_exchange, specialized, flavorbase):
#         if u in (mine-gus): 
#           is_not_there = False
#           break
#       if is_not_there: print "G", x, enumeration.as_bitstring(x, flavorbase)
#     for x in mine-gus:
#       is_not_there = True
#       for u, type in all_symmetrics(x, translations, label_exchange, specialized, flavorbase):
#         if u in (gus-mine): 
#           is_not_there = False
#           break
#       if is_not_there: print "M", x, enumeration.as_bitstring(x, flavorbase)
        nconf+=1 
t1 = time.time()
print "Took ", (t1 -t0)/60, "mn to complete."

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

  cell1 = atat.rMatrix3d( [ [  1.0,  0.5,  0.0 ], \
                            [  0.0,  0.5*sqrt(3.0),  0.0 ], \
                            [  0.0,  0.0,  1.0 ] ] )
  cell2 = atat.rMatrix3d( [ [  0.5,  0.5,  0.0 ], \
                            [  -0.5*sqrt(3.0),  0.5*sqrt(3.0),  0.0 ], \
                            [  0.0,  0.0,  1 ] ] )
  lattice.cell = cell1
  lattice.scale = 4.42

  lattice.sites.append( crystal.Site( (0., 0, 0) ) )
  lattice.sites[0].type = crystal.StringVector( [ "K", "Rb" ] );
  lattice.sites.append( crystal.Site( (0.5, sqrt(3)/6.0, 0.5) ) )
  lattice.sites[1].type = crystal.StringVector( [ "K", "Rb" ] );

  lattice.sites[0].pos = atat.rVector3d(0.5,  0.5/sqrt(3.0), 0.25)
  lattice.sites[1].pos = atat.rVector3d(0.5,  -0.5/sqrt(3.0), 0.75)

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

def enum( _n, _lattice ):
  from lada import enumeration, atat, crystal
  from math import pow
  import time

  supercells = enumeration.find_all_cells(lattice, _n)
  smiths = enumeration.create_smith_groups(_lattice, supercells)
  nflavors = enumeration.count_flavors(_lattice)
  nsites = len(_lattice.sites)
  transforms = enumeration.create_transforms(_lattice)
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
      cell = _lattice.cell * supercell.hermite
      specialized = []
      for transform in transforms:
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
        if specialized_database[x]: yield x, smith, supercell, flavorbase

# def main():
from lada import enumeration, atat, crystal
from math import pow
import time

lattice = create_lattice()
lattice.set_as_crystal_lattice()
species = ("K", "Rb")

nconf = 1
t0 = time.time()
nconf = 0
for n in range(2, 6):
  npern = 0
  oldsupercell = None
  for x, smith, supercell, flavorbase in enum(n, lattice):
    if oldsupercell == None or oldsupercell != supercell:
      npern = 0
    print "%5i %5i %2i " % (nconf, npern, n),
    for i in range(0,3):
      for j in range(0,i+1):
        print "%2i " % (supercell.hermite[i,j]),
    print "%10i %20s " % (int(x), enumeration.as_bitstring(x, flavorbase))
    npern += 1
    nconf += 1

t1 = time.time()
print "Took ", (t1 -t0)/60, "mn to complete."

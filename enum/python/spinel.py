#! /uhome/mdavezac//usr/bin/python
# 
#PBS -l nodes=1:ppn=1,walltime=00:30:0
#PBS -q Std
#PBS -m n 
#PBS -e err
#PBS -o out


def create_lattice():
  from lada import crystal, atat

  lattice = crystal.Lattice()

  lattice.cell = atat.rMatrix3d( [ [ 0, 0.5, 0.5 ], \
                                   [ 0.5, 0, 0.5 ], \
                                   [ 0.5, 0.5, 0 ] ] )

  # Manganese - Tetrahedral
  lattice.sites.append( crystal.Site( (0, 0, 0) ) )
  lattice.sites[0].type = crystal.StringVector( [ "Mg" ] );
  lattice.sites.append( lattice.sites[0] )
  lattice.sites[0].pos = atat.rVector3d( (   0,    0,    0) )
  lattice.sites[1].pos = atat.rVector3d( (0.25, 0.25, 0.25) )

  # Aluminum - Octahedral
  lattice.sites.append( crystal.Site( (5.0/8.0, 5.0/8.0, 5.0/8.0) ) )
  lattice.sites[2].type = crystal.StringVector( ["Al", "Mg"] )
  lattice.sites.extend( [ lattice.sites[2] for u in range(0,3) ] )
  lattice.sites[2].pos = atat.rVector3d( [5.0/8.0, 5.0/8.0, 5.0/8.0] )
  lattice.sites[3].pos = atat.rVector3d( [5.0/8.0, 7.0/8.0, 7.0/8.0] )
  lattice.sites[4].pos = atat.rVector3d( [7.0/8.0, 5.0/8.0, 7.0/8.0] )
  lattice.sites[5].pos = atat.rVector3d( [7.0/8.0, 7.0/8.0, 5.0/8.0] )

  # Oxygens
  x = 0.387
  lattice.sites.append( crystal.Site( (x, x, x) ) )
  lattice.sites[6].type = crystal.StringVector( ["O"] )
  lattice.sites.extend( [ lattice.sites[6] for u in range(0,7) ] )
  lattice.sites[ 6].pos = atat.rVector3d( [     x,      x,      x] )
  lattice.sites[ 7].pos = atat.rVector3d( [     x,     -x,     -x] )
  lattice.sites[ 8].pos = atat.rVector3d( [0.25-x, 0.25-x, 0.25-x] )
  lattice.sites[ 9].pos = atat.rVector3d( [0.25-x, 0.25+x, 0.25+x] )
  lattice.sites[10].pos = atat.rVector3d( [    -x,     -x,      x] )
  lattice.sites[11].pos = atat.rVector3d( [    -x,      x,     -x] )
  lattice.sites[12].pos = atat.rVector3d( [0.25+x, 0.25-x, 0.25+x] )
  lattice.sites[13].pos = atat.rVector3d( [0.25+x, 0.25+x, 0.25-x] )

# for site in lattice.sites: 
#   site.pos -= atat.rVector3d( (3.0/8.0, 3.0/8.0, 3.0/8.0) )

  lattice.scale = 7.5 # in bhor?
  lattice.find_space_group()

  return lattice

def link_sites(_lattice):
  from lada import crystal

  result = []
  accounted = []


  for i, site in enumerate(_lattice.sites):
    if i in accounted or len(site.type) == 1: continue
    
    result.append( [i,[]] )
    accounted.append(i)
    for op in _lattice.space_group:
      t = crystal.which_site( op(site.pos), _lattice )
      if t in accounted: continue

      accounted.append(t)
      result[-1][1].append( (t, op) )
  return result

def create_pairs(_lattice, _n):
  from lada import ce

  links = link_sites(_lattice)

  results = []
  for u in links:
    pairs = ce.create_pairs(_lattice, _n, u[0] )
    for p in pairs: print p
    results.append( (u[0], pairs) )
    for v in u[1]: 
      results.append( (v[0], pairs.apply_rotation(v[1])) )
  return results

def square_enum( _n, _lattice ):
  from lada import enumeration, atat, crystal
  from math import pow
  import numpy 
  import pyublas
  import time

  structure = crystal.sStructure()
  structure.cell = atat.rMatrix3d( [[_n, 0, 0], [0, _n, 0], [0, 0, _n]] )
  crystal.fill_structure(structure);
  smith = crystal.smith_normal_transform(structure) 

  nflavors = enumeration.count_flavors(_lattice)
  nsites = len([0 for u in _lattice.sites if  len(u.type) > 1])
  transforms = enumeration.create_transforms(_lattice)
  
  card = int(smith[1][0]*smith[1][1]*smith[1][2]*nsites)
  label_exchange=enumeration.LabelExchange( card, nflavors )
  flavorbase = enumeration.create_flavorbase(card, nflavors)
  translations = enumeration.Translation(smith[1], nsites)
  database = enumeration.Database(card, nflavors)
  maxterm = 0
  inv_cell = atat.inverse(_lattice.cell)
  for x in xrange(1, int(pow(nflavors, card))-1):
    if float( sum(enumeration.as_numpy(x, flavorbase)) ) != float(card) / 2.0:
      database[x] = False
      continue
      
    if not database[x]: continue
    maxterm = x
    for labelperm in label_exchange:
      t = labelperm(x, flavorbase)
      if t > x: database[t] = False
    for translation in translations:
      t = translation(x, flavorbase)
      if t > x: database[t] = False
      for labelperm in label_exchange:
        u = labelperm(t, flavorbase)
        if u > x: database[u] = False

    # checks supercell dependent transforms.
    mine = []
    # creates list of transformation which leave the supercell invariant.
    specialized = []
    for transform in transforms:
      if not transform.invariant(structure.cell): continue
      transform.init(smith[0], smith[1])
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
      if specialized_database[x]:
        yield x, inv_cell * structure.cell, flavorbase
        
def enum( _n, _lattice ):
  from lada import enumeration, atat, crystal
  from math import pow
  import numpy 
  import pyublas
  import time

  supercells = enumeration.find_all_cells(_lattice, _n)
  smiths = enumeration.create_smith_groups(_lattice, supercells)
  nflavors = enumeration.count_flavors(_lattice)
  nsites = len([0 for u in _lattice.sites if  len(u.type) > 1])
  transforms = enumeration.create_transforms(_lattice)
  for smith in smiths:
    card = smith.smith[0]*smith.smith[1]*smith.smith[2]*nsites
    label_exchange=enumeration.LabelExchange( card, nflavors )
    flavorbase = enumeration.create_flavorbase(card, nflavors)
    translations = enumeration.Translation(smith.smith, nsites)
    database = enumeration.Database(card, nflavors)
    maxterm = 0
    for x in xrange(1, int(pow(nflavors, card))-1):
      if float( sum(enumeration.as_numpy(x, flavorbase)) ) != float(card) / 2.0:
        database[x] = False
        continue
        
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
        if specialized_database[x]:
          yield x, supercell.hermite, flavorbase

def compute_sqs_parameters( _structure, _pairs ):
  from lada.ce import find_pis

  value = find_pis(_pairs[0][1], _structure, _pairs[0][0]); 
  for cls in _pairs[1:]:
    value += find_pis(cls[1], _structure, cls[0]); 
  return value
    

def main():
  from lada import crystal, ce, enumeration
  import random
  import sys
  import os
  import time

# sys.path.append( os.path.join(os.environ["HOME"], "usr/src/LaDa-3.30/enum/python/") )
  
# import enum

  # creates the lattice with smallest unit-cell.
  lattice = create_lattice()
  reduced_lattice = crystal.Lattice(lattice)
  reduced_lattice.sites.clear()
  for site in lattice.sites: 
    if len(site.type) > 1: reduced_lattice.sites.append(site)

  pairs = create_pairs(reduced_lattice, 2)
  lattice.set_as_crystal_lattice()

  nconf = 1
  t0 = time.time()
  nconf = 0
  values = []
  for n in range(2, 5):
    npern = 0
    oldsupercell = None
    structure = None
    for x, hermite, flavorbase in enum(n, lattice):
      if oldsupercell == None or oldsupercell != hermite:
        npern = 0
        structure = crystal.Structure()
        print structure.lattice()
        structure.cell = lattice.cell * hermite
        crystal.fill_structure(structure)
        print structure
        
#     print "%5i %5i %2i " % (nconf, npern, n),
#     for i in range(0,3):
#       for j in range(0,i+1):
#         print "%2i " % (supercell.hermite[i,j]),
#     print "%10i %20s " % (int(x), enumeration.as_bitstring(x, flavorbase))
      enumeration.as_structure(structure, x, flavorbase)
      dummy = compute_sqs_parameters(structure, pairs) 
      if dummy != None: print dummy
        # values.append( (dummy, hermite, x) )
      npern += 1
      nconf += 1

  return
  print "done computing sqs parameters."
  def compvals(_a, _b):
    from math import fabs
    a, b = 0, 0
    for u in _a[0]: a += fabs(u)
    for u in _b[0]: b += fabs(u)
    a /= float(len(_a[0]))
    b /= float(len(_a[0]))
    if fabs(a - b) < 1e-10: return 0
    if a < b: return -1
    return 1

  values = sorted(values, compvals) 
  for val in values[:10]:
    print val

  t1 = time.time()
  print "Took ", (t1 -t0)/60, "mn to complete."
  return


if __name__ == "__main__":
  main()

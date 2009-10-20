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
#     if float( sum(enumeration.as_numpy(x, flavorbase)) ) != float(card) / 2.0:
#       database[x] = False
        
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
          yield x, smith, supercell, flavorbase

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
  lattice.set_as_crystal_lattice()
  print lattice
  print link_sites(lattice)

  nconf = 1
  t0 = time.time()
  nconf = 0
  for n in range(2, 3):
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
  return


  # creates a superstructure with conventional unit-cell
  structure = crystal.sStructure()
  structure.cell = lattice.cell * atat.rMatrix3d( [ [ 1, 0, 0],
                                                    [ 0, 1, 0],
                                                    [ 0, 0, 1] ] )
  # fills in atoms in structure
  crystal.fill_structure( structure )

  # iterates over atoms in site 1 (Octahedral)
  for atom in filter(structure.atoms, lambda x: x.site == 1 ):
    # atom.site is the index of the site ad defined in the lattice.
    # atom.pos is the position.
    # atom.type is the type...
    atom.type = random.choice( ["Al", "Mg"] )
   
  
  # one way to print structure.
  print structure
  # prints in xcrysden format.
  print structure.xcrysden()
  # prints in POSCAR (VASP) format.
  crystal.print_poscar(structure)


if __name__ == "__main__":
  main()

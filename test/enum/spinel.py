
# look at inverse only.
def check_concentration(x, flavorbase, smith):
  from lada.enumeration import as_numpy

  types = as_numpy(x, flavorbase)
  return 3*len([i for i in types if i == 1]) == len(types)

def create_database(lattice, n0=1, n1=4):
  from lada.enumeration import Enum
 
  enum = Enum(lattice)
  # prints lattice
  print "lattice:\n  cell: ",
  for i in range(3):
    for j in range(3):
      print lattice.cell[i,j],
  print 
  print "  scale: ", lattice.scale
  for site in lattice.sites:
    print "  site: ",
    for i in range(3): print site.pos[i],
    for i in site.type: print i, 
    print 
  print "endlattice\n\n"

  nconf = 1
  nconf = 0
  result = []
  perrow = 21
  iter = 0
  for n in range(n0, n1):
    if iter != 0: 
      print 
      iter = 0
    print "\nn: ", n
    npern = 0
    oldsupercell = None
    first = True
    for x, smith, supercell, flavorbase in enum.xn(n, check_concentration):
      if oldsupercell == None or oldsupercell != supercell:
        if first:
          print "  flavorbase: ", len(flavorbase), flavorbase[1]
          first = False
        npern = 0
        if iter != 0: 
          print 
          iter = 1
        print "  hermite: ",
        for i in range(3):
          for j in range(3):
            print supercell.hermite[i,j],
        print
        print "  transform: ",
        for i in range(3):
          for j in range(3):
            print supercell.transform[i,j],
        print
        print "  smith: ",
        for i in range(3): print smith.smith[i],
        print
        oldsupercell = supercell
        iter = 0

      if iter == 0:
        print "   ",
      elif iter % perrow == 0:
        iter = 0
        print
        print "   ",
      iter += 1
      print x,

def main():
  from lada.crystal.A2BX4 import b5

  # create spinel lattice.
  lattice = b5()
  for site in lattice.sites:
    if "X" in site.type: continue
    site.type = 'A', 'B'

  create_database(lattice, 1, 5)

main()


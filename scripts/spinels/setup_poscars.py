
#! /uhome/mdavezac/usr/bin/python
# 
#PBS -l nodes=2:ppn=2,walltime=48:00:0
#PBS -q Std
#PBS -m n 
#PBS -e err
#PBS -o out

def create_structure_dir(structure, prefix=None):
  from os import makedirs
  from os.path import join, exists

  # create structure name
  nAl = len([ 0 for v in structure.atoms if v.type == "Al" ])
  nMg = len([ 0 for v in structure.atoms if v.type == "Mg" ])
  nO = len([ 0 for v in structure.atoms if v.type == "O" ])
  structure.name = "Al%iMg%iO%i" % ( nAl, nMg, nO )

  # create directory
  dirname = structure.name
  if prefix != None: dirname = join(prefix, dirname)
  u = 1
  while exists(dirname) :
    dirname = "%s-%i" % (structure.name, u) 
    if prefix != None: dirname = join(prefix, dirname)
    u += 1

  # create directory, POSCAR, and sendme 
  makedirs(dirname)
  return dirname

def choose_structures( howmany, filename="database", prefix=None ):
  from random import randint
  from lada import crystal, enumeration
  import database

  N = 0
  for structure in database.read_database(filename, withperms=True): N  += 1

  which = set([])
  while len(which) != howmany: which.add(randint(0,N-1))
  which = sorted(list(which))
  for n, structure in enumerate(database.read_database(filename, withperms=True)):
    if n != which[0]: continue

    which.pop(0)

    dirname = create_structure_dir(structure, prefix)
    structure.name += ", database %s, structure %i" % (filename, n)
    crystal.print_poscar(structure, ("Al", "Mg", "O"), dirname )

    if len(which) == 0: break
    

if __name__ == "__main__":
  choose_structures(35)

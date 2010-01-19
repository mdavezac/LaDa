
#! /uhome/mdavezac/usr/bin/python
# 
#PBS -l nodes=2:ppn=2,walltime=48:00:0
#PBS -q Std
#PBS -m n 
#PBS -e err
#PBS -o out
import os

def create_structure_dir(structure):
    # create structure name
    nAl = len([ 0 for v in structure.atoms if v.type == "Al" ])
    nMg = len([ 0 for v in structure.atoms if v.type == "Mg" ])
    nO = len([ 0 for v in structure.atoms if v.type == "O" ])
    structure.name = "Al%iMg%iO%i" % ( nAl, nMg, nO )

    # create directory
    dirname = structure.name
    u = 1
    while os.path.exists(dirname) :
      dirname = "%s-%i" % (structure.name, u) 
      u += 1

    # create directory, POSCAR, and sendme 
    os.mkdir(dirname)
    return dirname

def choose_structures( howmany, filename="database" ):
  import os
  import shutil
  import random
  from lada import crystal, enumeration
  import database

  N = 0
  for structure in database.read_database(filename, withperms=True): N  += 1

  which = set([])
  while len(which) != howmany:
    which.add(random.randint(0,N-1))
  which = sorted(list(which))
  for n, structure in enumerate(database.read_database(filename, withperms=True)):
    if n != which[0]: continue

    which.pop(0)

    dirname = create_structure_dir(structure)
    structure.name += ", database %s, structure %i" % (filename, n)
    crystal.print_poscar(structure, ("Al", "Mg", "O"), dirname )

    if len(which) == 0: break
    

if __name__ == "__main__":
  choose_structures(35)

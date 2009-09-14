#! /usr/bin/python

def extract( _dir, _vasp ):
  from os import getcwd, walk, path
  from vasp import Vasp
  
  pwd  = getcwd()
  for root, dirs, files in walk(_dir):
    _vasp.indir = path.join(_dir, root)
    try: yield _vasp.final()
    except GeneratorExit : return
    except: pass

def input_directories(_filename):

  import re
  import sys
  from os import path

  file = 0
  file = open(_filename, "r")

  for line in file: 
    if len(line) == 0: continue
    if re.search( r"^\s*\#", line ): continue
    for dir in line.split():
      if dir == "#": break
      if not path.exists(dir):
        print dir, "does not exist."
        continue
      if not path.isdir(dir):
        print dir, "is not a valid directory."
        continue
      yield dir 

def main():
  def create_vasp():
    from vasp import Vasp, Specie
    K = Specie( "K", "~/AtomicPotentials/pseudos/K_s" )
    Rb = Specie( "Rb", "~/AtomicPotentials/pseudos/Rb_s" )

    vasp = Vasp()
    vasp.species = [K, Rb]
    vasp.cutoff_safety
    vasp.encut  = 1.3 * vasp.__get_encut__( vasp.species )
    return vasp

  import re
  import sys
  from os import path

  if "-h" in sys.argv or "--help" in sys.argv:
    print ">", sys.argv[0], " filename [default: input]" 
    print "Extracts output from directories listed in filename."
    print "Each line of in filename is a directory containing other directories "\
          "with static vasp calculations."

  filename = "input"
  if len(sys.argv) == 2: filename = str(sys.argv[1])
  elif len(sys.argv) > 2: 
    print "Too many parameters on command-line."
    sys.exit(1)

  for dir in input_directories(filename):
    for structure in extract(dir, vasp):
      print structure


if __name__ == "__main__":
  main();
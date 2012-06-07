def test(program):
  from numpy import all, abs
  from tempfile import mkdtemp
  from shutil import rmtree
  from gafunc import Functional as GAFunc
  from lada import default_comm
  # something screwing with the seed... 
  # possibly uuid. Can't make this a deterministic test so easily.
  functional = GAFunc(program, 10, popsize=20, rate=0.5)
  dir = '/tmp/test' #mkdtemp()
  try: 
    result = functional(dir, default_comm)
    indiv =  result.best(1)
    assert all(abs(indiv.genes-1) < 1e-8)
  finally:
    pass

   #try: rmtree(dir)
   #except: pass

if __name__ == "__main__":
  from sys import argv, path
  from os.path import abspath
  if len(argv) < 1: raise ValueError("test need to be passed location of pifunc.")
  if len(argv) > 2: path.extend(argv[2:])
  test(abspath(argv[1]))

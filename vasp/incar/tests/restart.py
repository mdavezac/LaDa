def test():
  from pickle import loads, dumps
  from os import chdir, getcwd, remove
  from os.path import join, exists
  from shutil import rmtree
  from tempfile import mkdtemp
  from lada.vasp.incar._params import Restart
  from lada.vasp.files import WAVECAR, CHGCAR
  from restart_class import Extract, Vasp

  texts = "ISTART = 0   # start from scratch\n"\
          "ICHARG = 2   # superpositions of atomic densities", \
          "ISTART = 0   # start from scratch\n"\
          "ICHARG = 12   # superpositions of atomic densities", \
          "ISTART = 0   # start from scratch\n"\
          "ICHARG = 2   # superpositions of atomic densities"
  directory = mkdtemp()
  try: 
    for e in [None, Extract(directory, False)]:
      for v, t in zip([None, Vasp(True), Vasp(False)], texts):
        assert Restart(e).incar_string(vasp=v) == t
        r = loads(dumps(Restart(e)))
        assert r.incar_string(vasp=v) == t
        if hasattr(e, 'directory'):
          assert r.value.directory == directory 
          assert r.value.success == e.success
    
    with open(join(directory, CHGCAR), "w") as file: pass
    with open(join(directory, WAVECAR), "w") as file: pass
    for e in [None, Extract(directory, True)]:
      for v, t in zip([None, Vasp(True), Vasp(False)], texts):
        assert Restart(e).incar_string(vasp=v) == t
        r = loads(dumps(Restart(e)))
        assert r.incar_string(vasp=v) == t
        if hasattr(e, 'directory'):
          assert r.value.directory == directory 
          assert r.value.success == e.success
    
    cwd = getcwd()
    texts = "ISTART = 1   # restart\n"\
            "ICHARG = 1   # from charge {0}".format(join(directory, CHGCAR)),\
            "ISTART = 1   # restart\n"\
            "ICHARG = 11   # from charge {0}".format(join(directory, CHGCAR)),\
            "ISTART = 1   # restart\n"\
            "ICHARG = 1   # from charge {0}".format(join(directory, CHGCAR))
    with open(join(directory, CHGCAR), "w") as file: file.write('hello')
    for e in [Extract(directory, True)]:
      for v, t in zip([None, Vasp(True), Vasp(False)], texts):
        playdir = mkdtemp()
        try: 
          chdir(playdir)
          assert Restart(e).incar_string(vasp=v) == t
          r = loads(dumps(Restart(e)))
          assert r.incar_string(vasp=v) == t
          if hasattr(e, 'directory'):
            assert r.value.directory == directory 
            assert r.value.success == e.success
          assert exists(join(playdir, CHGCAR))
        finally:
          chdir(cwd)
          rmtree(playdir)
  
    texts = "ISTART = 1   # restart\n"\
            "ICHARG = 0   # from wavefunctions {0}".format(join(directory, WAVECAR)),\
            "ISTART = 1   # restart\n"\
            "ICHARG = 10   # from wavefunctions {0}".format(join(directory, WAVECAR)),\
            "ISTART = 1   # restart\n"\
            "ICHARG = 0   # from wavefunctions {0}".format(join(directory, WAVECAR))
    with open(join(directory, WAVECAR), "w") as file: file.write('hello')
    for e in [Extract(directory, True)]:
      for v, t in zip([None, Vasp(True), Vasp(False)], texts):
        playdir = mkdtemp()
        try: 
          chdir(playdir)
          assert Restart(e).incar_string(vasp=v) == t
          r = loads(dumps(Restart(e)))
          assert r.incar_string(vasp=v) == t
          if hasattr(e, 'directory'):
            assert r.value.directory == directory 
            assert r.value.success == e.success
          assert exists(join(playdir, WAVECAR))
          assert exists(join(playdir, CHGCAR))
        finally:
          chdir(cwd)
          rmtree(playdir)
  
    remove(join(directory, CHGCAR))
    with open(join(directory, CHGCAR), "w") as file: pass
    for e in [Extract(directory, True)]:
      for v, t in zip([None, Vasp(True), Vasp(False)], texts):
        playdir = mkdtemp()
        try: 
          chdir(playdir)
          assert Restart(e).incar_string(vasp=v) == t
          r = loads(dumps(Restart(e)))
          assert r.incar_string(vasp=v) == t
          if hasattr(e, 'directory'):
            assert r.value.directory == directory 
            assert r.value.success == e.success
          assert exists(join(playdir, WAVECAR))
          assert not exists(join(playdir, CHGCAR))
        finally:
          chdir(cwd)
          rmtree(playdir)
  
  finally: rmtree(directory)

if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()


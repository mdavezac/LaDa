def test(executable):
  """ Tests JobFolderProcess. Includes failure modes.  """
  from tempfile import mkdtemp
  from os.path import join, exists
  from shutil import rmtree
  from numpy import all, arange, abs, array
  from lada.jobfolder.jobfolder import JobFolder
  from lada.jobfolder.massextract import MassExtract
  from lada.jobfolder import save
  from lada.process.jobfolder import JobFolderProcess
  from lada.process import Fail
  from lada.misc import Changedir
  from lada import default_comm
  from functional import Functional

  root = JobFolder()
  for n in xrange(8):
    job = root / str(n)
    job.functional = Functional(executable, [n])
    job.params['sleep'] = 1

  comm = default_comm.copy()
  comm['n'] = 4

  dir = mkdtemp()
  try: 
    program = JobFolderProcess(root, nbpools=2, outdir=dir, comm=comm)
    assert program.nbjobsleft > 0
    program.start()
    assert len(program.process) == 2
    while not program.poll():  continue
    assert program.nbjobsleft == 0
    save(root, join(dir, 'dict.dict'), overwrite=True)
    extract = MassExtract(join(dir, 'dict.dict'))
    assert all(extract.success.itervalues())
    order = array(extract.order.values()).flatten()
    assert all(arange(8) - order == 0)
    pi = array(extract.pi.values()).flatten()
    assert all(abs(pi - array([0.0, 3.2, 3.162353, 3.150849,
                               3.146801, 3.144926, 3.143907, 3.143293]))\
                < 1e-5 )
    error = array(extract.error.values()).flatten()
    assert all(abs(error - array([3.141593, 0.05840735, 0.02076029, 0.009256556,
                                  0.005207865, 0.00333321, 0.002314774, 0.001700664]))\
                < 1e-5 )
    assert all(n['n'] == comm['n'] for n in extract.comm)
    # restart
    assert program.poll()
    assert len(program.process) == 0
    program.start()
    assert len(program.process) == 0
    assert program.poll()
  finally:
    try: rmtree(dir)
    except: pass

  try: 
    job = root / str(666)
    job.functional = Functional(executable, [666])
    program = JobFolderProcess(root, nbpools=2, outdir=dir, comm=comm)
    assert program.nbjobsleft > 0
    program.wait()
    assert program.nbjobsleft == 0
  except Fail: 
    assert len(program.errors.keys()) == 1
    assert '666' in program.errors
  else: raise Exception
  finally:
    try: rmtree(dir)
    except: pass
  try: 
    job.functional.order = [667]
    program = JobFolderProcess(root, nbpools=2, outdir=dir, comm=comm)
    assert program.nbjobsleft > 0
    program.wait()
    assert program.nbjobsleft == 0
  finally:
    try: rmtree(dir)
    except: pass

if __name__ == "__main__":
  from sys import argv, path
  from os.path import abspath
  if len(argv) < 1: raise ValueError("test need to be passed location of pifunc.")
  if len(argv) > 2: path.extend(argv[2:])
  test(abspath(argv[1]))

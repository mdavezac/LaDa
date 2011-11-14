#! python
""" Runs one job from the jobtree. """
def main():
  import re 
  from sys import path as python_path
  from os import getcwd, environ
  from os.path import abspath, join, relpath, exists
  from argparse import ArgumentParser
  from lada import jobs
  from lada.opt.changedir import Changedir

  # below would go additional imports.

  parser = ArgumentParser( prog="runone", description = re.sub("\\s+", " ", __doc__[1:]))
  parser.add_argument( "--jobid", dest="n", default=1, \
                       help="Job number", metavar="N", type=int )
  parser.add_argument( "--relative", dest="relative", default=None, \
                       help="Perform calculations in a directory relative "
                            "current, but starting at RELATIVE, rather than HOME.",
                       metavar="RELATIVE" )
  parser.add_argument( "--ppath", dest="ppath", default=None, \
                       help="Directory to add to python path",
                       metavar="Directory" )
  parser.add_argument('--external', action="store_true", dest="external", \
                      help="Launches jobs as external program, not library. Only for VASP at this point.")
  parser.add_argument('--nprocs', dest="nprocs", default=1, type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--timeout', dest="timeout", default=300, type=int,\
                      help="Time to wait for job-dictionary to becom available "
                           "before timing out (in seconds). A negative or null "
                           "value implies forever. Defaults to 5mn.")
  parser.add_argument('pickle', metavar='FILE', type=str, help='Path to a jobdictionary.')

  try: options = parser.parse_args()
  except SystemExit: return

  # is workdir relative
  if options.relative is not None: 
    # get path relative to home.
    if options.relative not in environ:
      print "Error: could not find environment variable", options.relative, "."
      print "Will work on default dir instead."
      options.relative = None

  # additional path to look into.
  if options.ppath is not None: python_path.append(options.ppath)

  if not exists(options.pickle): 
    print "Could not find file {0}.".format(options.pickle)
    return

  # whether to launch programs externally or internally.
  if options.external: 
    from lada.mpi import NullComm
    comm = NullComm(options.nprocs)
  else:
    from lada.mpi import world
    comm = world
  timeout = None if options.timeout <= 0 else options.timeout
  
  # loop over all jobs -- Needs communicator!
  jobdict = jobs.load(options.pickle, comm=comm)
  for i, (outdir, job) in enumerate(jobdict.iteritems()):
    if i != options.n: continue
    if options.relative is None: 
      job.compute(comm=comm, outdir=outdir, inplace=True, external=options.external)
    else:
      # Computes relative path... made complicated by cray's compute vs head node setup.
      workdir = abspath(outdir)
      with Changedir(environ["HOME"]) as cwd:
        workdir = join(environ[options.relative], relpath(workdir, getcwd()))
      # now pass on relative workdir, where HOME is substituted by options.relative.
      job.compute(comm=comm, outdir=outdir, workdir=workdir, inplace=False, external=options.external)

if __name__ == "__main__":
  from sys import argv
  if "--external" in argv:
    import lada
    lada.lada_with_mpi = False
    main()
  else:
    import traceback
    from boost.mpi import Exception as mpiException, abort
    from lada.mpi import world
    try: main()
    except mpiException as e: 
      if world.is_root: 
        traceback.print_exc()
        file.stderr("Encountered mpi exception %s."% (e))
        print 
      abort(0)
      raise
    except: # other exceptions
      if world.is_root: 
        traceback.print_exc()
        print 
      abort(0)
      raise

#! python
""" Runs one job from the jobtree. """
def main():
  import cPickle
  import re 
  from sys import path as python_path
  from os import getcwd, environ
  from os.path import expanduser, abspath, join, relpath
  from optparse import OptionParser
  from boost.mpi import world
  from lada import jobs
  from lada.opt.changedir import Changedir

  # below would go additional imports.

  parser = OptionParser( description = re.sub("\\s+", " ", __doc__[1:]),\
                         usage = "usage: %prog [options] filename_1 .. filename_n" )
  parser.add_option( "--jobid", dest="n", default=1, \
                     help="Job number", metavar="N", type="int" )
  parser.add_option( "--relative", dest="relative", default=None, \
                     help="Perform calculations in a directory relative "
                          "current, but starting at RELATIVE, rather than HOME.",
                     metavar="RELATIVE" )
  parser.add_option( "--ppath", dest="ppath", default=None, \
                     help="Directory to add to python path",
                     metavar="Directory" )

  try: (options, args) = parser.parse_args()
  except SystemExit: return

  # is workdir relative
  if options.relative != None: 
    # get path relative to home.
    if options.relative not in environ:
      print "Error: could not find environment variable", options.relative, "."
      print "Will work on default dir instead."
      options.relative = None

  # additional path to look into.
  if options.ppath != None: python_path.append(options.ppath)

  if len(args) == 0:
    print "No pickle specified on input. Eg, need a filename on input."
    return

  # loop over all jobs -- Needs communicator!
  jobdict = jobs.load(args[0],comm=world)
  for i, (job, outdir) in enumerate(jobdict.iteritems()):
    if i != options.n: continue
    if options.relative == None: 
      out = job.compute(comm=world, outdir=outdir, inplace=True)
    else:
      # Computes relative path... made complicated by cray's compute vs head node setup.
      workdir = abspath(outdir)
      with Changedir(environ["HOME"]) as cwd:
        workdir = join(environ[options.relative], relpath(workdir, getcwd()))
      # now pass on relative workdir, where HOME is substituted by options.relative.
      out = job.compute(comm=world, outdir=outdir, workdir=workdir, inplace=False)

if __name__ == "__main__":
  from sys import stderr
  import traceback
  from boost.mpi import Exception as mpiException, abort, world
  try: main()
  except mpiException as e: 
    for i in range(world.size): 
      if i == world.rank:
        traceback.print_exc()
        file.stderr("Encountered mpi exception %s."% (e))
        print 
      world.barrier()
    abort(0)
    raise
  except: # other exceptions
    for i in range(world.size): 
      if i == world.rank:
        traceback.print_exc()
        print 
      world.barrier()
    abort(0)
    raise


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

  (options, args) = parser.parse_args()

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

  # loop over all jobs
  print "getcwd ", getcwd()
  for i, (job, outdir) in enumerate(jobs.walk_through(args[0])):
    if i != options.n: continue
    if options.relative == None: 
      print "outdir: ", outdir
      out = job.compute(comm=world, outdir=outdir, inplace=True)
    else:
      # Computes relative path... made complicated by cray's compute vs head node setup.
      workdir = abspath(outdir)
      with Changedir(environ["HOME"]) as cwd:
        workdir = join(environ[options.relative], relpath(workdir, getcwd()))
      # now pass on relative workdir, where HOME is substituted by options.relative.
      print "outdir: ", outdir
      print "workdir: ", workdir
      out = job.compute(comm=world, outdir=outdir, workdir=workdir, inplace=False)

if __name__ == "__main__":
  from sys import stderr
  import traceback
  from boost.mpi import Exception as mpiException, abort
  try: main()
  except mpiException as e: 
    file.stderr("Encountered mpi exception %s."% (e))
    traceback.print_exc()
    abort(0)
  except: # other exceptions
    traceback.print_exc()
    abort(0)
    raise


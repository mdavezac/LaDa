#! python
""" Runs a job-tree, with optional parallelization. 

    The jobtree is read from the file given on input. It should be a pickled
    JobDict instance. If more than one file is given, then the jobtree is
    concatenation of each jobtree pickled in each file. Parallelization can
    occur over pbs-scripts, in which case calculations executed in other
    pbs-scripts are ignored here (as given by npbs and pbspools command line
    arguments).  Parallelization can also occur over mpi-processes.
""" 
def main():
  import cPickle
  import re 
  from optparse import OptionParser
  from boost.mpi import world
  from lada import jobs

  # below would go additional imports.

  parser = OptionParser( description = re.sub("\\s+", " ", __doc__[1:]),\
                         usage = "usage: %prog [options] filename_1 .. filename_n" )
  parser.add_option( "--pbspools", dest="pbspools", help="Number of pbspools",\
                     metavar="N", type="int", default=1 )
  parser.add_option( "--npbs", dest="npbs", help="Which pbs script is this", metavar="N",\
                     type="int", default=0)
  parser.add_option( "--procpools", dest="procpools", default=1, \
                     help="Number of process pools per job", metavar="N", type="int" )

  # below would go additional program options.

  (options, args) = parser.parse_args()
  assert options.pbspools > options.npbs

  if isinstance(args, str): jobpickle = [args]
  elif len(args) == 0:
    print "No pickle specified on input. Eg, need a filename on input."
    return
  elif len(args) == 1: jobpickle = args

  # below would go additional out-of-loop code.
    
  # create jobs.
  jobs.current = jobs.JobDict()
  for file in jobpickle: jobs.current.update( jobs.load(path=file) )

  # makes sure that the number of pools is not too large.
  if world.size < options.procpools: options.procpools = world.size
  # creates local comms.
  color = world.rank %  options.procpools
  local_comm = world.split(color)

  # total numbe of separate pools.
  totpools = options.pbspools * options.procpools
  # index in total number of pools
  n = options.npbs * options.procpools + color
  # loop over all jobs
  for i, (job, outdir) in enumerate(jobs.walk_through()):
    # bypasses those jobs not done here.
    if i % totpools != n: continue
    out = job.compute(comm = local_comm, outdir = outdir)
    # below would go additional inner loop code.
if __name__ == "__main__": main()

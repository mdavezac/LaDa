#! python
""" Runs one job from the jobtree. """
def main():
  import re 
  from sys import path as python_path
  from copy import deepcopy
  from os.path import exists
  from argparse import ArgumentParser
  from lada import jobs, default_comm

  # below would go additional imports.

  parser = ArgumentParser( prog="runone", description = re.sub("\\s+", " ", __doc__[1:]))
  parser.add_argument( "--jobid", dest="names", default='+', \
                       help="Job name", metavar="N", type=str )
  parser.add_argument( "--ppath", dest="ppath", default=None, \
                       help="Directory to add to python path",
                       metavar="Directory" )
  parser.add_argument('--nbprocs', dest="nbprocs", default=default_comm['n'], type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--ppn', dest="ppn", default=default_comm['ppn'], type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--timeout', dest="timeout", default=300, type=int,\
                      help="Time to wait for job-dictionary to becom available "
                           "before timing out (in seconds). A negative or null "
                           "value implies forever. Defaults to 5mn.")
  parser.add_argument('pickle', metavar='FILE', type=str, help='Path to a job-folder.')

  try: options = parser.parse_args()
  except SystemExit: return

  # additional path to look into.
  if options.ppath is not None: python_path.append(options.ppath)

  if not exists(options.pickle): 
    print "Could not find file {0}.".format(options.pickle)
    return

  # whether to launch programs externally or internally.
  comm = deepcopy(default_comm)
  comm['n'] = options.nbprocs
  comm['ppn'] = options.ppn
  timeout = None if options.timeout <= 0 else options.timeout
  
  jobfolder = jobs.load(options.pickle, timeout=timeout)
  for name in options.names:
    jobfolder[name].compute(comm=comm, outdir=name[1:], external=options.external)

if __name__ == "__main__": main()

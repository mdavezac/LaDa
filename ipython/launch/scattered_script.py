""" Runs one job from the jobfolder. """
def main():
  import re 
  from sys import path as python_path
  from os.path import exists
  from argparse import ArgumentParser
  from pylada import jobfolder
  from pylada.process.mpi import create_global_comm
  import pylada

  # below would go additional imports.

  parser = ArgumentParser( prog="runone", description = re.sub("\\s+", " ", __doc__[1:]))
  parser.add_argument( "--jobid", dest="names", nargs='+', type=str, \
                       help="Job name", metavar="N" )
  parser.add_argument( "--ppath", dest="ppath", default=None, \
                       help="Directory to add to python path",
                       metavar="Directory" )
  parser.add_argument('--nbprocs', dest="nbprocs", default=pylada.default_comm['n'], type=int,\
                      help="Number of processors with which to launch job.")
  parser.add_argument('--ppn', dest="ppn", default=pylada.default_comm['ppn'], type=int,\
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

  # Set up mpi processes.
  pylada.default_comm['ppn'] = options.ppn
  pylada.default_comm['n'] = options.nbprocs
  create_global_comm(options.nbprocs)

  timeout = None if options.timeout <= 0 else options.timeout
  
  jobfolder = jobfolder.load(options.pickle, timeout=timeout)
  for name in options.names:
    jobfolder[name].compute(comm=pylada.default_comm, outdir=name)

if __name__ == "__main__": main()

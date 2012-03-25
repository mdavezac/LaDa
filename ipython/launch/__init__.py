""" IPython launch magic function. """
__docformat__ = "restructuredtext en"

def launch(self, event):
  """ Launches PBS/slurm jobs.


      The usual setup is to first explore a dictionary of some sort, 
      then modify it and save it, and finally launch it.

      >>> # to recompute errors.
      >>> %explore errors path/to/original/pickle
      >>> %goto next
      >>> # modify dictionary here.
      >>> %showme functional
      >>> ...
      >>> %goto next
      >>> %showme functional
      >>> ...
      >>> # Saves the modified job.
      >>> # the new path could a different filename in the same directory.
      >>> # This way, the unsuccessful output from the first run will be
      >>> # overwritten.
      >>> %savejobs path/to/original/pickle.errors
      >>> # then  launch.
      >>> %launch scattered --walltime "24:00:00"
  """ 
  import argparse
  from ...jobs import load as load_jobs
  from ... import interactive, default_comm
# from .scattered import parser as scattered_parser
  from .interactive import parser as interactive_parser

  # main parser
  parser = argparse.ArgumentParser(prog='%launch')
  # options supported by all.
  opalls = argparse.ArgumentParser(add_help=False)
  opalls.add_argument( 'pickle', metavar='FILE', type=str, nargs='*', default="", 
                       help='Optional path to a jobdictionary. If not present, the '\
                            'currently loaded job-dictionary will be launched.')
  opalls.add_argument( '--nbprocs', type=int, default=default_comm.get('n', 1),
                       nargs='?', help="Number of processes over which to launch calculations." )
  opalls.add_argument( '--force', action="store_true", dest="force", \
                       help="Launches all untagged jobs, even those "\
                            "which completed successfully." )

  # subparsers
  subparsers = parser.add_subparsers(help='Launches one job per untagged calculations')

  # launch scattered.
# scattered_parser(self, subparsers, opalls) 
  interactive_parser(self, subparsers, opalls) 

  # parse arguments
  try: args = parser.parse_args(event.split())
  except SystemExit as e: return None
  
  comm = default_comm.copy()
  comm['n'] = args.nbprocs

  # creates list of dictionaries.
  jobdicts = []
  if args.pickle != '':
    for pickle in args.pickle:
      try: d = load_jobs(path=pickle, timeout=20)
      except ImportError as e:
        print "ImportError: ", e
        return
      except Exception as e:
        print e
        if LockFile(pickle).is_locked:
          print "You may want to check for the existence of {0}."\
                .format(LockFile(pickle).lock_directory)
          print "If you are sure there are no jobs out there accessing {0},\n"\
                "you may want to delete that directory.".format(args.pickle)
          return
      else: jobdicts.append((d, pickle))
  else: # current job dictionary.
    if interactive.jobdict is None:
      print "No current job-dictionary."
      return
    if interactive.jobdict_path is None:
      print "No path for currrent job-dictionary."
      return
    jobdicts = [(interactive.jobdict, interactive.jobdict_path)]
  
  # calls specialized function.
  args.func(self, args, jobdicts, comm)



def completer(self, info):
  """ Completion for launchers. """
  from .scattered import completer as scattered_completer
  from .interactive import completer as interactive_completer
  from IPython import TryNext
  from .. import jobdict_file_completer

  data = info.line.split()[1:]
  if "scattered" in data: return scattered_completer(self, info, data)
  elif "interactive" in data: return interactive_completer(self, info, data)
  return ["scattered", "interactive"]
         

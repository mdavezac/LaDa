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
  from ...jobfolder import load as load_jobs
  from ... import interactive, default_comm
  from .scattered import parser as scattered_parser
  from .interactive import parser as interactive_parser
  from ...misc import RelativePath

  # main parser
  parser = argparse.ArgumentParser(prog='%launch')
  # options supported by all.
  opalls = argparse.ArgumentParser(add_help=False)
  opalls.add_argument( 'pickle', metavar='FILE', type=str, nargs='*', default="", 
                       help='Optional path to a job-folder. If not present, the '\
                            'currently loaded job-dictionary will be launched.')
  opalls.add_argument( '--force', action="store_true", dest="force", \
                       help="If present, launches all untagged jobs, even those "\
                            "which completed successfully." )
  # subparsers
  subparsers = parser.add_subparsers(help='Launches one job per untagged calculations')

  # launch scattered.
  scattered_parser(self, subparsers, opalls) 
  interactive_parser(self, subparsers, opalls) 

  # parse arguments
  try: args = parser.parse_args(event.split())
  except SystemExit as e: return None
  
  # creates list of dictionaries.
  jobfolders = []
  if args.pickle != '':
    for pickle in args.pickle:
      pickle = RelativePath(pickle).path
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
      else: jobfolders.append((d, pickle))
  else: # current job folder.
    if interactive.jobfolder is None:
      print "No current job-dictionary."
      return
    if interactive.jobfolder_path is None:
      print "No path for currrent job-dictionary."
      return
    jobfolders = [(interactive.jobfolder, interactive.jobfolder_path)]
  
  # calls specialized function.
  args.func(self, args, jobfolders)



def completer(self, info):
  """ Completion for launchers. """
  from .scattered import completer as scattered_completer
  from .interactive import completer as interactive_completer
  from IPython import TryNext
  from .. import jobfolder_file_completer

  data = info.line.split()[1:]
  if "scattered" in data: return scattered_completer(self, info, data)
  elif "interactive" in data: return interactive_completer(self, info, data)
  return ["scattered", "interactive"]
         

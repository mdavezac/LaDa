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
  from .. import _get_current_job_params
  from .scattered import parser as scattered_parser
  from .interactive import parser as interactive_parser

  # main parser
  parser = argparse.ArgumentParser(prog='%launch')
  # options supported by all.
  opalls = argparse.ArgumentParser(add_help=False)
  opalls.add_argument( 'pickle', metavar='FILE', type=str, nargs='*', default="", 
                       help='Optional path to a jobdictionary. If not present, the '\
                            'currently loaded job-dictionary will be launched.')
  opalls.add_argument( '--external', action="store_true", dest="external", \
                       help="Launches jobs as external program, not library." )


  # subparsers
  subparsers = parser.add_subparsers(help='Launches one job per untagged calculations')

  # launch scattered.
  scattered_parser(self, subparsers, opalls) 
  interactive_parser(self, subparsers, opalls) 

  # parse arguments
  try: args = parser.parse_args(event.split())
  except SystemExit as e: return None

  # creates list of dictionaries.
  pickles = set(args.pickle) - set([""])
  if len(pickles) > 0: 
    jobdicts = []
    for p in pickles:
      try: d = load_jobs(path=p)
      except Exception as e: 
        print "JobDict could not be loaded form {0}.\n{e}\n".format(p,e=e)
        return
      jobdicts.append((d, p))
  else: # current job dictionary.
    current, path = _get_current_job_params(self, 2)
    if current == None: return
    if path == None: return
    jobdicts = [(current, path)]
  
  # calls specialized function.
  args.func(self, args, jobdicts)



def completer(self, info):
  """ Completion for launchers. """
  from .scattered import completer as scattered_completer
  from .interactive import completer as interactive_completer
  from IPython.ipapi import TryNext

  data = info.line.split()
  types = ["scattered", "interactive"]
  if len(data)  <= 2 and data[-1] not in types: return types 
  if data[1] == "scattered": return scattered_completer(self, info, data)
  if data[1] == "interactive": return interactive_completer(self, info, data)
  raise TryNext
         

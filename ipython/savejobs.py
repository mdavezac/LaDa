""" Save a job to a dictionary. """
def savejobs(self, event):
  """ Saves current job to current filename and directory. """
  from os.path import exists, abspath, isfile
  from ..jobs import JobParams, MassExtract as Collect, save as savejobs
  from lada import interactive
  jobdict = interactive.jobdict.root
  jobdict_path = interactive.jobdict_path

  if jobdict is None:
    print "No job-dictionary to save." 
    return
  args = [u for u in event.split() ]
  if len(args) > 1: 
    print "savejobs takes one or no arguments."
    return

  if len(args) == 1: jobdict_path = args[0]

  if jobdict_path is None: 
    print "No current job-dictionary path.\n"\
          "Please specify on input, eg\n"\
          ">saveto this/path/filename"
    return
  if exists(jobdict_path): 
    if not isfile(jobdict_path): 
      print "{0} is not a file.".format(jobdict_path)
      return
    a = ''
    while a not in ['n', 'y']:
      a = raw_input("File {0} already exists.\nOverwrite? [y/n] ".format(jobdict_path))
    if a == 'n':
      print "Aborting."
      return
  savejobs(jobdict.root, jobdict_path, overwrite=True, timeout=10) 
  if len(args) == 1:
    if "collect" not in self.user_ns: self.user_ns["collect"] = Collect(dynamic=True, path=jobdict_path)
    if "jobparams" not in self.user_ns: self.user_ns["jobparams"] = JobParams()

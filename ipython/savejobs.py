""" Save a job to a dictionary. """
def savejobs(self, event):
  """ Saves current job to current filename and directory. """
  from os.path import exists, abspath, isfile
  from ..jobs import JobParams, MassExtract as Collect, save as savejobs
  from .. import interactive
  from ..misc import RelativePath
  jobfolder = interactive.jobfolder.root
  jobfolder_path = interactive.jobfolder_path

  if jobfolder is None:
    print "No job-dictionary to save." 
    return
  args = [u for u in event.split() ]
  if len(args) > 1: 
    print "savejobs takes one or no arguments."
    return

  if len(args) == 1:
    jobfolder_path = RelativePath(args[0]).path
    interactive.jobfolder_path = jobfolder_path

  if jobfolder_path is None: 
    print "No current job-folder path.\n"\
          "Please specify on input, eg\n"\
          ">saveto this/path/filename"
    return
  if exists(jobfolder_path): 
    if not isfile(jobfolder_path): 
      print "{0} is not a file.".format(jobfolder_path)
      return
    a = ''
    while a not in ['n', 'y']:
      a = raw_input("File {0} already exists.\nOverwrite? [y/n] ".format(jobfolder_path))
    if a == 'n':
      print "Aborting."
      return
  savejobs(jobfolder.root, jobfolder_path, overwrite=True, timeout=10) 
  if len(args) == 1:
    if "collect" not in self.user_ns: self.user_ns["collect"] = Collect(dynamic=True, path=jobfolder_path)
    if "jobparams" not in self.user_ns: self.user_ns["jobparams"] = JobParams()

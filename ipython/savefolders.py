def savefolders(self, event):
  """ Saves job-folder to disk.
  
      This function can be called in one of three ways:

      >>> savefolders filename.dict rootfolder 
  
      In this case, "filename.dict" is a file where to save the jobfolder
      "rootfolder". The latter must be a python variable, not another filename. 
      The current job-folder becomes "rootfolder", and the current path 


      >>> savefolders filename.dict

      Saves the current job-folder to "filename.dict". Fails if no current
      job-folder.

      >>> savefolders 

      Saves the current job-folder to the current job-folder path. Fails if
      either are unknown.
  """
  from os.path import exists, isfile
  from ..jobfolder import JobParams, MassExtract as Collect, save
  from .. import interactive
  from ..misc import RelativePath

  args = [u for u in event.split()]
  if '--help' in args or '-h' in args:
    print savefolders.__doc__
    return 

  if len(args) > 2: 
    print "savefolders takes zero, one, or two arguments."
    return

  if len(args) == 2:
    from .explore import explore
    explore(self, args[1])
    savefolders(self, args[0])
    return

  if interactive.jobfolder is None: 
    print "No current job-folder."
    print "Please load first with %explore."
    return
  jobfolder = interactive.jobfolder.root
  jobfolder_path = interactive.jobfolder_path

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
  save(jobfolder.root, jobfolder_path, overwrite=True, timeout=10) 
  if len(args) == 1:
    if "collect" not in self.user_ns: self.user_ns["collect"] = Collect(dynamic=True, path=jobfolder_path)
    if "jobparams" not in self.user_ns: self.user_ns["jobparams"] = JobParams()

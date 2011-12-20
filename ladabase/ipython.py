""" IPython magic functions to interact with the OUTCAR database. """

""" List of files previously pushed to ladabase. """
def get_file_list(self, args):
  """ Returns list of files to act upon. """
  from os.path import isfile, join, isdir
  from glob import iglob
  from itertools import chain
  from ..vasp import MassExtract
  files = []
  for var in chain(*(iglob(v) for v in args.outcar)):
    if isfile(var): # case where this is a file.
      files.append(var)
    elif isdir(var): # case where this is a directory.
      for extract in MassExtract(var).values():
        files.append(join(extract.directory, extract.OUTCAR))
    else: raise ValueError("Could not make sense of {0}.".format(var))
  return files


def push(self, cmdl):
  """ Pushes directory with OUTCARs or single OUTCAR to the database. """
  try: from .. import ladabase_root_push
  except: return 
  try: from .. import username as pymongo_username
  
  except:
    print "Could not find username. Please edit the file '~/.lada', and add:"
    print ">>> username = \"Jane Doe\""
    print "Without '>>>', with 'username' flushed left, and your name "\
          "within explicit quotation marks on the right hand side." 
    return 

  
  import argparse
  import tarfile 
  from datetime import datetime
  from os import getcwd
  from os.path import relpath, join
  from getpass import getusername
  from IPython.ipapi import TryNext
  from ..vasp import Extract
  from .misc import get_username, get_ladabase
  import re

  try: ladabase = get_ladabase()
  except RuntimeError as e: print e; return; 
  try: get_username()
  except RuntimeError as e: print e; return; 
  
  parser = argparse.ArgumentParser(prog='%push',
                     description='Push single OUTCAR or directory of OUTCARS to ladabase. ')
  parser.add_argument( 'outcar', metavar='OUTCAR', type=str, nargs='*',
                       help='Job dictionary, OUTCAR file, or root directory.' )
  parser.add_argument( 'showdir', action="store_true",
                       help="Print directory name where OUTCAR archives "\
                            " and comment file are stored." )
  try: args = parser.parse_args(cmdl.split())
  except SystemExit as e: return None

  if args.showdir: 
    print ladabase_root_push
    return 

  # list all structures.
  files, errors = [], ""
  for file in get_file_list(self, args):
    if file in ladabase: print "File {0} is already in the database.".format(file)
    extract = Extract(file)
    if not extract.success:
      errors += "**** File {0} is not a successful calculation.\n".format(file)
    else: files.append((file, extract))
  if len(files) == 0: 
    print "Nothing to add to database."
    return

  # Now creates an archive.
  filename = getusername() + str(datetime.now()) 
  filename = join(ladabase_root_push, filename)
  tarme = tarfile.open(filename + ".tgz", 'w:gz')
  directory = getcwd()

  # gets a comment to go with the push.
  notefile = filename + ".comment"
  with open(notefile, "w") as file:
    file.write('\n# files simultaneously pushed to the database:\n')
    for filename, extract in files: file.write('#   ' + relpath(filename, getcwd()) + '\n')
  try: self.shell.hooks.editor(notefile, 0)
  except TryNext:
    print "Could not open editor."
    return
  with open(notefile, "r") as file: comment = file.read()
  stripped = re.sub('#.*\n', '', comment)
  if len(stripped.replace('\n','').replace(' ', '')) == 0: 
    print "Empty comment. Aborting."
    return
  # Now adds stuff add end of comment.
  with open(notefile, 'w') as file:
    file.write(stripped + "\n")
    for path, extract in file:
      file.write("# file: {0}\n".format(relpath(path, directory)))
    file.write("# operator: {0}\n".format(pymongo_username))

  for file, extract in files:
    tarme.add(file, arcname=relpath(file, directory))
    print "Pushing ", relpath(file, directory)
  tarme.close()
  print 
  print errors

#   with open(file, 'r') as outcar:
#     kwargs = {'compression': compression, 'comment':comment}
#     kwargs['is_dft'] =  extract.is_dft
#     kwargs['is_gw'] =  extract.is_gw
#     added = ladabase.push( relpath(file, getcwd()), outcar.read(), **kwargs)
#     if added is not None:
#       just_added.append(added)
#       print "Pushed {0}.".format(file, getcwd())

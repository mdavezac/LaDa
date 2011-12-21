""" IPython magic functions to interact with the OUTCAR database. """

def _getcomment(self, filename = None):
  """ Gets comment from user. """
  from sys import stderr
  from os import fdopen
  from socket import gethostname
  from IPython.ipapi import TryNext
  try: from .. import username
  except ImportError:
    print >>stderr, "Could not import username with which to tag files in database.\n"\
                    "Please add `username = 'my full name'` in ~/.lada.\n"
    return
  if len(username) == 0:
    print >>stderr, "Username with which to tag files in database is empty.\n"\
                    "Please add `username = 'my full name'` in ~/.lada.\n"
    return

  import re
  if filename is None:
    from tempfile import mkstemp
    # gets a comment to go with the push.
    notefile, filename = mkstemp()
    fdopen(notefile).close()
  with open(filename, 'w') as file:
    file.write("\n# operator: {0}\n# hostname: {1}".format(username, gethostname()))
  try: self.shell.hooks.editor(filename, 0)
  except TryNext:
    print "Could not open editor."
    return
  with open(filename, "r") as file: comment = file.read()
  stripped = re.sub('#.*\n', '', comment)
  if len(stripped.replace('\n','').replace(' ', '')) == 0: 
    print "Empty comment. Aborting."
    return
  with open(filename, 'a') as file:
    file.write("\n# operator: {0}\n# hostname: {1}".format(username, gethostname()))
  return filename

def _walk_files(args):
  """ Loops over valid OUTCAR files. """
  from os import walk
  from os.path import isfile, isdir, join
  from glob import iglob
  from re import compile, search
  from ..vasp import Extract

  excludes = [compile('relax_cellshape'), compile('relax_ions')]
  if args.exclude is not None:
    for i in args.exclude: excludes.append(compile(i))

  for input in args.directories:
    if any(search(i, input) is not None for i in excludes): continue
    for path in iglob(input):
      if any(search(i, path) is not None for i in excludes): continue
      if isfile(path): 
        extract = Extract(path)
        if not extract.success: continue
        if extract.relaxation != "static": continue
        yield path
      elif isdir(path): 
        for root, dirs, files in walk(path):
          if args.norecurrence: dirs[:] = []
          else: dirs[:] = [d for d in dirs if all(search(i, d) is None for i in excludes)]
          if args.pattern is not None: 
            for file in iglob(join(root, args.pattern)):
              if any(search(i, file) is not None for i in excludes): continue
              extract = Extract(file)
              if not extract.success: continue
              if extract.relaxation != "static": continue
              yield file
          else:
            for file in files:
              file = join(root, file)
              if any(search(i, file) is not None for i in excludes): continue
              extract = Extract(file)
              if not extract.success: continue
              if extract.relaxation != "static": continue
              yield file

def _get_local_push_parser():
  """ Returns parser for local push. """
  import argparse
  parser = argparse.ArgumentParser(prog='%push',
                     description='Push single OUTCAR or directory of OUTCARS to ladabase. ')
  parser.add_argument( 'directories', metavar='DIR/OUTCAR', type=str, nargs='*',
                       help='Job dictionary, OUTCAR file, or root directory.' )
  parser.add_argument( '--list', action="store_true", help='List files which will be pushed and exits.' )
  parser.add_argument( '--norecurrence', action="store_true", help='Do not walk into subdirectory.')
  parser.add_argument( '--pattern', type=str, nargs=1, default=None,
                       help='When searching within directories, examine only files '\
                            'fitting this pattern. Defaults to examining all files.')
  parser.add_argument( '--exclude', type=str, nargs='*', help="Files and directories matching "\
                       "these regex patterns are excluded from the search.")
  return parser

def local_push(self, cmdl):
  """ Pushes file to a local directory using tar.
  
      This process is daemonized.
  """
  from os import makedirs
  from os.path import join, exists, dirname
  from getpass import getuser
  from datetime import datetime
  from subprocess import call as syscall, STDOUT
  from sys import executable
  from .. import  local_push_dir


  try: args = _get_local_push_parser().parse_args(cmdl.split())
  except SystemExit: return None
  if args.list:
    for u in _walk_files(args): print u
    return 

  if not exists(local_push_dir): makedirs(local_push_dir)
  filename = join(local_push_dir, getuser() + str(datetime.now()).replace(' ', '_'))

  # gets comment file.
  comment_filename = _getcomment(self, filename + ".comment")
  if comment_filename is None: return

  commands = ["nohup", executable, join(dirname(__file__), "_local_tar.py"),
              '--tarfile', filename, '--commentfile', comment_filename ]
  commands.extend(cmdl.split())
  syscall(commands, stdout=open(filename + ".log", 'w'), stderr=STDOUT)

  print "Working in background. "
  print "See", (filename + ".log"), "for files being pushed, or you use --list option."
  print "You can log out at this point."

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


def push_database(self, cmdl):
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
  from getpass import getuser
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
  filename = getuser() + str(datetime.now()) 
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
    with open(file, 'r') as outcar:
      kwargs = {'compression': 'bz2', 'comment':comment}
      kwargs['is_dft'] =  extract.is_dft
      kwargs['is_gw'] =  extract.is_gw
      if ladabase.push( relpath(file, getcwd()), outcar.read(), **kwargs) is not None:
        print "Pushed {0}.".format(file, getcwd())
  tarme.close()
  print 
  print errors


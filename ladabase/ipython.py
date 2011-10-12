""" IPython magic functions to interact with the OUTCAR database. """

""" List of files previously pushed to ladabase. """

def get_file_list(self, args):
  """ Returns list of files to act upon. """
  from os.path import exists, isfile, join, isdir
  from ..vasp import MassExtract
  files = []
  for var in args.outcar:
    if exists(var): # case where this is a directory or file.
      if isfile(var): # case where this is a file.
        files.append(var)
      elif isdir(var): # case where this is a directory.
        for extract in MassExtract(var).values():
          files.append(join(extract.directory, extract.OUTCAR))
      else: raise ValueError("Could not make sense of {0}.".format(var))
    elif var in self.api.user_ns: 
      for extract in self.api.user_ns[var]: files.append(join(extract.directory, extract.OUTCAR))
  return files



def push(self, cmdl):
  """ Pushes directory with OUTCARs or single OUTCAR to the database. """
  import argparse
  from os import getcwd
  from os.path import relpath
  from IPython.ipapi import TryNext
  from ..vasp import Extract
  from ..record import Record
  import re

  if 'ladabase' not in self.api.user_ns:
    print "Could not find ladabase instance."
    return
  ladabase = self.api.user_ns['ladabase']
  
  parser = argparse.ArgumentParser(prog='%push',
                     description='Push single OUTCAR or directory of OUTCARS to ladabase. ')
  parser.add_argument( 'outcar', metavar='OUTCAR', type=str, default="current_jobdict", nargs='*',
                       help='Job dictionary, OUTCAR file, or root directory.' )
  parser.add_argument( '--compression', type=str, default="none", dest="compression", nargs='?',
                       help='Type of compression used.' )

  try: args = parser.parse_args(cmdl.split())
  except SystemExit as e: return None

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


  # gets a comment to go with the push.
  notefile = self.shell.mktempfile()
  with open(notefile, "w") as file:
    file.write('\n# files simultaneously pushed to the database:\n')
    for filename, extract in files: file.write('#   ' + relpath(filename, getcwd()) + '\n')

  try: self.shell.hooks.editor(notefile, 0)
  except TryNext:
    print "Could not open editor."
    return
  with open(notefile, "r") as file: comment = file.read()
  stripped = re.sub('#.*\n', '\n', comment)
  if len(stripped.replace('\n','').replace(' ', '')) == 0: 
    print "Empty comment. Aborting."
    return

  compression = None if args.compression.lower() == "none" else args.compression

  just_added = []
  for file, extract in files: 
    with open(file, 'r') as outcar:
      kwargs = {'compression': compression, 'comment':comment}
      kwargs['is_dft'] =  extract.is_dft
      kwargs['is_gw'] =  extract.is_gw
      added = ladabase.push( relpath(file, getcwd()), outcar.read(), **kwargs)
      if added is not None:
        just_added.append(added)
        print "Pushed {0}.".format(file, getcwd())

  if len(just_added) == 0: return
  record = Record()
  if hasattr(record, '_pushed'):
    _pushed = record._pushed
    _pushed.append((just_added, comment, compression))
    record._pushed = _pushed
  else: record._pushed = [(just_added, comment, compression)]
  if len(errors) != 0: print "\n", errors


def amend(self, cmdl):
  """ Pushes directory with OUTCARs or single OUTCAR to the database. """
  from IPython.ipapi import TryNext
  from hashlib import sha512
  from re import sub as replace
  from ..record import Record

  if 'ladabase' not in self.api.user_ns:
    print "Could not find ladabase instance."
    return
  ladabase = self.api.user_ns['ladabase']
  
  if len(cmdl.replace(' ', '')) != 0:
    print "amend does not require command-line arguments."
    return
  # list all structures.
  record = Record()
  if not hasattr(record, '_pushed'): 
    print 'Nothing was pushed from this directory in living memory.'
    return
  ids, comment, compression = record._pushed[-1]
  hash = sha512(comment.replace(' ', '').replace('\n', '')).hexdigest()

  notefile = self.shell.mktempfile()
  with open(notefile, "w") as file: file.write(comment)
  try: self.shell.hooks.editor(notefile, 0)
  except TryNext:
    print "Could not open editor."
    return
  with open(notefile, "r") as file: comment = file.read()
  if hash == sha512(comment.replace(' ', '').replace('\n', '')).hexdigest():
    print "No change made to comment. Aborting."
    return
  stripped = replace('#.*\n', '\n', comment)
  if len(stripped.replace('\n','').replace(' ', '')) == 0: 
    print "Empty comment. Aborting."
    return

  for id in ids: 
    if len(list(ladabase.files.find({'_id': id}))) == 0: 
      print "Could not find object with id {0}".format(str(id))
    ladabase.files.update({'_id': id}, {'$set': {'comment': comment}})
  _pushed = record._pushed
  _pushed[-1] = (ids, comment, compression)
  record._pushed = _pushed

  





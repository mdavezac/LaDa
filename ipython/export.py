""" IPython export magic function. """
__docformat__ = "restructuredtext en"


def export(self, event):
  """ Tars files from a calculation.  """
  import argparse
  import tarfile
  from os import getcwd
  from os.path import exists, isfile, extsep, relpath
  from ..opt import RelativeDirectory

  parser = argparse.ArgumentParser(prog='%export',
                     description='Exports input/output files from current jobdictionary. '\
                                 'Depending on the extension of FILE, this will create '\
                                 'a simple tar file, or a compressed tar file. Using the '\
                                 'option --list, one can also obtain a list of all files '\
                                 'which would go in the tar archive. '\
                                 'Finally, this function only requires the \"collect\" '\
                                 'exists in the usernamespace. It may have been declared '\
                                 'from loading a jobdictionary using \"explore\", or directly '\
                                 'with \"collect = vasp.MassExtract()\".' )
  parser.add_argument( 'filename', metavar='FILE', type=str, default='export.tar.gz',
                       nargs='?',
                       help='Path to the tarfile. Suffixes ".gz" and ".tgz" indicate '\
                            'gzip compression, whereas ".bz" and ".bz2" indicate bzip '\
                            'compression. Otherwise, no compression is used.')
  parser.add_argument( '--incar', action="store_true", dest=" incar",
                       help='Include INCAR files.' )
  parser.add_argument( '--doscar', action="store_true", dest="doscar",
                       help='Include DOSCAR files.' )
  parser.add_argument( '--poscar', action="store_true", dest="poscar",
                       help='Include POSCAR files.' )
  parser.add_argument( '--chgcar', action="store_true", dest="chgcar",
                       help='Include CHGCAR files.' )
  parser.add_argument( '--contcar', action="store_true", dest="contcar",
                       help='Include POTCAR files.' )
  parser.add_argument( '--potcar', action="store_true", dest="potcar",
                       help='Include POTCAR files.' )
  parser.add_argument( '--wavecar', action="store_true", dest="wavecar",
                       help='Include WAVECAR files.' )
  parser.add_argument( '--list', action="store_true", dest="aslist",
                       help='Do not tar, return a list of all the files.' )

  try: args = parser.parse_args(event.split())
  except SystemExit as e: return None

  if 'collect' not in self.api.user_ns:
    print "Could not find 'collect' object in user namespace."
    print "Please load a job-dictionary."
    return

  kwargs = args.__dict__.copy()
  kwargs.pop('filename', None)
  jobdict = self.api.user_ns.get('current_jobdict_path', None)
  if args.aslist:
    from IPython.genutils import SList
    try: 
      result = SList([file for file in self.api.user_ns['collect'].iterfiles(**kwargs)])
      if jobdict != None: result.append(jobdict)
      return result
    except Exception as e:
      print "Encoutered error while tarring file."
      print e
  else:
    args.filename = relpath(RelativeDirectory(args.filename).path, getcwd())
    if exists(args.filename): 
      if not isfile(args.filename):
        print "{0} exists but is not a file. Aborting.".format(args.filename)
        return 
      a = ''
      while a not in ['n', 'y']:
        a = raw_input("File {0} already exists.\nOverwrite? [y/n] ".format(args.filename))
      if a == 'n': print "Aborted."; return

    if args.filename.find(extsep) == -1: endname = ''
    else: endname = args.filename[-args.filename[::-1].find(extsep)-1:][1:]
    if endname in ['gz', 'tgz']:   tarme = tarfile.open(args.filename, 'w:gz')
    elif endname in ['bz', 'bz2']: tarme = tarfile.open(args.filename, 'w:bz2')
    else:                          tarme = tarfile.open(args.filename, 'w')
    try: 
      for file in self.api.user_ns['collect'].iterfiles(**kwargs): tarme.add(file)
      if jobdict != None: tarme.add(jobdict)
    except Exception as e:
      print "Encoutered error while tarring file."
      print e
    tarme.close()
    print "Saved archive to {0}.".format(args.filename)

def completer(self, event):
  """ Completer for export. """
  from glob import iglob
  from itertools import chain
  from os.path import isdir
  data = set(event.line.split()) - set(["export", "%export"])
  result = set(['--incar', '--doscar', '--poscar', '--chgcar', '--contcar', 
                '--potcar', '--wavecar', '--list'])
  other = data - result
  if len(other) == 0: 
    for file in chain( iglob('*.tar'), iglob('*.tar.gz'), 
                       iglob('*.tgz'), iglob('*.bz'), iglob('*.bz2') ):
      result.add(file)
    result |= set([u for u in self.api.magic("%mglob dir:*")])
  elif len(other) == 1: 
    other = other.__iter__().next()
    string = "{0}*.tar {0}*.tar.gz {0}*.tgz {0}*.tar.bz {0}*.tar.bz2 dir:{0}*".format(other)
    result |= set([u for u in self.api.magic("%mglob " + string)])
    if isdir(other) and other[-1] != '/':
      string = "{0}/*.tar {0}/*.tar.gz {0}/*.tgz {0}/*.tar.bz "\
               "{0}/*.tar.bz2 dir:{0}/*".format(other)
      result |= set([u for u in self.api.magic("%mglob " + string)])
  return list(result - data)



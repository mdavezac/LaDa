""" Actually performs local transfer 

    As a stand-alone script so we can call from Popen with nohup.
"""
from sys import stderr, exit
from os import fork, setsid, getcwd, remove
from os.path import relpath
from shutil import move
from lada.ladabase.ipython import _walk_files, _get_local_push_parser
import tarfile 

# first fork to orphan the process created in the second.
parser = _get_local_push_parser()
parser.add_argument('--tarfile', type=str)
parser.add_argument('--commentfile', type=str)

try: args = parser.parse_args()
except SystemExit: exit(0)


try: 
  pid = fork()
  if pid > 0: exit(0)
except OSError as error: 
  print >>stderr, "Fork failed when daemonizing local push: {0.errno} ({0.strerror}.".format(error)
  exit(1)

# let first forked process grab session. (I think).
setsid()

# second fork to be orphaned.
try: 
  pid = fork()
  if pid > 0: exit(0)
except OSError as error: 
  print >>stderr, "Fork failed when daemonizing local push: {0.errno} ({0.strerror}.".format(error)
  exit(1)

# now gets list of files.
tarme = tarfile.open(args.tarfile + "_partial" + ".tgz", 'w:gz')
tarme.add(args.commentfile, arcname="comment")
cwd = getcwd()
for file in _walk_files(args):
  print file
  tarme.add(file, arcname=relpath(file, cwd))
tarme.close()

move(args.tarfile + "_partial" + ".tgz", args.tarfile + ".tgz")
# remove(args.commentfile)

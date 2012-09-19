""" Provides a quick launch of jmol. """
def jmol(self, event):
  """ Launches JMOL on a structure. """
  from inspect import ismethod
  from tempfile import NamedTemporaryFile
  from subprocess import call
  from os import remove
  from numpy import all, abs, max, min
  from ..crystal import write
  from .. import jmol_program
  from . import get_shell

  if '-h' in event.split() or '--help' in event.split(): 
    print  "usage: %jmol [-h] structure\n"                                     \
           "\n"                                                                \
           "Launches jmol for a given structure.\n"                            \
           "\n"                                                                \
           "positional arguments:\n"                                           \
           "  structure   Variable/expression referencing a structure.\n"      \
           "\n"                                                                \
           "optional arguments:\n"                                             \
           "  -h, --help  show this help message and exit"
    return 
  if len(event.rstrip().lstrip()) == 0:
    print '%jmol requires at least one argument.'
    return
  shell = get_shell(self)
  if event.rstrip().lstrip() in shell.user_ns:
    structure = shell.user_ns[event.rstrip().lstrip()]
  else: 
    structure = eval(event.rstrip().lstrip(), shell.user_ns.copy())
  if ismethod(getattr(structure, 'eval', None)):
    structure = structure.eval()

  if all(abs(structure.cell[:, 2] - [0, 0, 500.0]) < 1e-8):
    mini = abs(min([a.pos[2] for a in structure]))
    maxi = abs(max([a.pos[2] for a in structure]))
    structure.cell[2, 2] = mini + maxi
    structure.cell[2, 2] *= 1.1

  if all(abs(structure.cell[:, 1] - [0, 500.0, 0.0]) < 1e-8):
    mini = abs(min([a.pos[1] for a in structure]))
    maxi = abs(max([a.pos[1] for a in structure]))
    structure.cell[1, 1] = mini + maxi
    structure.cell[1, 1] *= 1.1

  if all(abs(structure.cell[:, 0] - [500.0, 0.0,  0.0]) < 1e-8):
    mini = abs(min([a.pos[0] for a in structure]))
    maxi = abs(max([a.pos[0] for a in structure]))
    structure.cell[0, 0] = mini + maxi
    structure.cell[0, 0] *= 1.1

  with NamedTemporaryFile(delete=False) as file:
    name = file.name
    file.write(write.castep(structure))
  try: call(jmol_program + ' ' + name, shell=True)
  except KeyboardInterrupt: pass
  finally: remove(name)

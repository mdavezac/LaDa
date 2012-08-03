""" Python wrapper for the CRYSTAL dft code. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Extract', 'AffineTransform', 'DisplaceAtoms', 'InsertAtoms',
            'Marker', 'ModifySymmetry', 'RemoveAtoms', 'Slabcut', 'read',
            'Slabinfo','Crystal', 'Molecule', 'Slab', 'Functional', 'Shell',
            'read_gaussian_basisset', 'MassExtract' ]

from .basis import Shell
from .functional import Functional
from .extract import Extract, MassExtract
from .geometry import AffineTransform, DisplaceAtoms, InsertAtoms, Marker,      \
                      ModifySymmetry, RemoveAtoms, Slabcut, Slabinfo, Slab
from .crystal import Crystal
from .molecule import Molecule
registered = { 'atomrot':    AffineTransform,
               'atomdisp':   DisplaceAtoms,
               'atominse':   InsertAtoms,
               'marker':     Marker,
               'modisymm':   ModifySymmetry,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'crystal':    Crystal,
               'molecule':   Molecule }
""" Keywords used in creating the geometry. """

def read(path):
  """ Reads CRYSTAL input from file

      Reads a file and finds the CRYSTAL_ input within. It then parses the file
      and creates a :py:class:`~functional.Functional` and
      :py:class:`~molecule.Molecule` (or derived) object.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      :returns: A two-tuple containing the structure and the functional

      :rtype: (:py:class:`~molecule.Molecule`, :py:class:`~functional.Functional`)

      .. _CRYSTAL: http://www.crystal.unito.it/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from .parse import parse
  b = parse(path)
  key = b.keys()[0]

  func = Functional()
  crys = Crystal()
  func.read_input(b)
  crys.read_input(b[key]['CRYSTAL'])
  crys.name = key
  return crys, func


def read_gaussian_basisset(path):
  """ Reads a GAUSSIAN94 basis set

      This is meant to read input files from the `basis set exchange website`__
      in the GAUSSIAN94_ format.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      .. __: https://bse.pnl.gov/bse/portal
      .. _GAUSSIAN94: http://www.gaussian.com/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from ..error import GrepError
  from ..misc import RelativePath
  if isinstance(path, str): 
    if path.find('\n') == -1:
      with open(RelativePath(path).path) as file: return read_gaussian_basisset(file)
    else:
      return read_gaussian_basisset(path.split('\n').__iter__())

  for line in path: 
    if set(['*']) == set(line.rstrip().lstrip()): break

  # read specie type
  try: line = path.next().split()
  except StopIteration: raise GrepError('Unexpected end of file')
  specie = line[0]
  result = {specie: []}
  
  try: 
    while True:
      line = path.next().split()
      if len(line) != 3:
        line = path.next().split()
        if len(line) != 2: break
        specie = line[0]
        result[specie] = []
        continue
      type, n, scale = line[0].lower(), int(line[1]), float(line[2])
      shell = Shell(type)
      for i in xrange(n): shell.append(*path.next().split())
      result[specie].append(shell)

  except StopIteration: pass

  return result

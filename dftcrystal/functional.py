__docformat__ = "restructuredtext en"
from .input import AttrBlock
from ..functools import stateless, assign_attributes
from .extract import Extract

class Functional(AttrBlock):
  """ Wrapper for the CRYSTAL program. """
  Extract = Extract
  """ Extraction class. """
  def __init__(self, copy=None, program=None, **kwargs):
    """ Creates the crystal wrapper. """
    from .hamiltonian import Dft
    from .basis import BasisSet
    from .optgeom import OptGeom
    from .input import TypedKeyword
    super(Functional, self).__init__()
    self.basis   = BasisSet()
    """ Holds definition of basis functions. """
    self.dft     = Dft()
    """ Holds definition of functional. """
    self.optgeom = OptGeom()
    """ Holds definition of geometry optimization. """
    self.title   = None
    """ Title of the calculation. 
    
        Overriden by the name of the input structure, if it exists.
    """
    self.program = program
    """ Path to crystal program.

        If this attribute is None, then :py:data:`~lada.crystal_program` is
        used.
    """ 
    self.add_keyword('maxcycle', TypedKeyword('maxcycle', int))
    self.add_keyword('tolinteg', TypedKeyword('tolinteg', [int]*6))
    self.add_keyword('toldep', TypedKeyword('toldep', int))
    self.add_keyword('tolpseud', TypedKeyword('tolpseud', int))
    self.add_keyword('toldee', TypedKeyword('toldee', int))

  def read_input(self, tree):
    """ Reads file or string with CRYSTAL input. """
    from ..error import IOError
    from .. import CRYSTAL_geom_blocks as starters

    self.title = tree.keys()[0]
    tree = tree[self.title]
    # read optgeom bit.
    found = False
    for starter in starters:
      if starter in tree.keys(): found = True; break
    if found == False:
      raise IOError('Could not find start of input in file.')
    if 'OPTGEOM' in tree[starter].keys():
      self.optgeom.read_input(tree[starter]['OPTGEOM'])

    # read basis set
    if 'BASISSET' in tree.keys(): 
      self.basis.read_input(tree['BASISSET'])

    # read hamiltonian stuff.
    if 'HAMSCF' in tree.keys(): 
      others = AttrBlock()
      others.read_input(tree['HAMSCF'])
      dft = others._crysinput.pop('DFT', None)
      if dft is not None: self.dft.read_input(dft)
      self._crysinput.update(others._crysinput)

  def print_input(self, **kwargs):
    """ Dumps CRYSTAL input to string. """

    # create optgeom part first, since it needs be inserted in the structure
    # bit. Split along lines and remove empty lines at end.
    # if empty, then make it empty.
    optgeom = self.optgeom.print_input(**kwargs).rstrip().split('\n')
    while len(optgeom[-1].rstrip().lstrip()) == 0: optgeom.pop(-1)
    if len(optgeom) == 2: optgeom = []

    result = ''
    if 'structure' in kwargs:
      structure = kwargs['structure']
      # insert name of structure as title.
      if hasattr(structure, 'name'):
        result += str(structure.name).rstrip().lstrip() + '\n'
      elif getattr(self, 'title', None) is not None:
        result += str(self.title).rstrip().lstrip()
      else: result += '\n'
      result = result.rstrip()
      if result[-1] != '\n': result += '\n'
       
      # figure out input of structure. remove empty lines.
      lines = structure.print_input(**kwargs).split('\n')
      while len(lines[-1].rstrip().lstrip()) == 0: lines.pop(-1)
      # insert optgeom inside the geometry bit.
      if len(optgeom) > 0: lines.insert(-2, optgeom)
      # turn back into string.
      result += '\n'.join(lines)
    else: # no structures. Not a meaningful input, but whatever.
      result += '\n'.join(optgeom)

    # now add basis
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    result += self.basis.print_input(**kwargs)

    # now add hamiltonian
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    result += self.dft.print_input(**kwargs)

    # add keywords contained directly in functional.
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    a = AttrBlock()
    a._crysinput = self._crysinput
    result += a.print_input(**kwargs)

    # end input and return
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    return result + 'END\n'

  def guess_workdir(outdir):
    """ Tries and guess working directory. """
    from os import environ, getpid
    from os.path import join
    from datetime import datetime
    return join( environ.get('PBS_TMPDIR', outdir),
                 '{0!s}.{1}'.format(datetime.today(), getpid()))

  def bringup(self, structure, workdir, restart):
    """ Creates file environment for run. """
    from os.path import join
    from ..misc import copyfile, Changedir
    from .. import CRYSTAL_filenames

    with Changedir(workdir) as cwd:
      # first copies file from restart
      if restart is not None: 
        for key, value in CRYSTAL_filenames.iteritems():
          copyfile( value.format('crystal'), key, nocopyempty=True,
                    symlink=False, nothrow="never" )
      # then copy files from current directory.
      for key, value in CRYSTAL_filenames.iteritems():
        copyfile( join(restart.directory, value.format('crystal')), 
                  key, nocopyempty=True, symlink=False, 
                  nothrow="never" )

      # then creates input file.
      string = self.print_input(crystal=self, structure=structure)
      with open('crystal.d12', 'w') as file: file.write(string)

  def bringdown(self, structure, workdir, outdir):
    """ Copies files back to output directory. 
    
        Prefixes output with input.
    """
    from os.path import join
    from ..misc import copyfile, Changedir
    from .. import CRYSTAL_filenames

    with Changedir(outdir) as cwd:
      for key, value in CRYSTAL_filenames.iteritems():
        copyfile( join(workdir, key), value.format('crystal'),
                  nocopyempty=True, symlink=False, nothrow="never" )

      with open('crystal.d12', 'r') as file: input = file.read()
      with open('crystal.out', 'r') as file: output = file.read()
      header = ''.join(['#']*20)
      with open('crystal.out', 'w') as file:
        file.write('{0} {1} {0}\n'.format(header, 'INPUT FILE'))
        input = input.rstrip()
        if input[-1] != '\n': input += '\n'
        file.write(input)
        file.write('{0} END {1} {0}\n'.format(header, 'INPUT FILE'))
        file.write(output)
        file.write('\n{0} {1} {0}\n'.format(header, 'FUNCTIONAL'))
        file.write(repr(self))
        file.write('\n{0} END {1} {0}\n'.format(header, 'FUNCTIONAL'))
        file.write('\n{0} {1} {0}\n'.format(header, 'STRUCTURE'))
        file.write(repr(structure))
        file.write('\n{0} END {1} {0}\n'.format(header, 'STRUCTURE'))

  
  @assign_attributes(ignore=['overwrite', 'comm', 'workdir'])
  @stateless
  def iter(self, structure, outdir=None, workdir=None, comm=None,
           overwrite=False, **kwargs):
    """ Performs a vasp calculation 
     
        If successfull results (see :py:attr:`extract.Extract.success`) already
        exist in outdir, calculations are not repeated. Instead, an extraction
        object for the stored results are given.

        :param structure:  
            :py:class:`~lada.crystal.Structure` structure to compute, *unless*
            a CONTCAR already exists in ``outdir``, in which case this
            parameter is ignored. (This feature can be disabled with the
            keyword/attribute ``restart_from_contcar=False``).
        :param outdir:
            Output directory where the results should be stored.  This
            directory will be checked for restart status, eg whether
            calculations already exist. If None, then results are stored in
            current working directory.
        :param comm:
            Holds arguments for executing CRYSTAL externally.
        :param overwrite:
            If True, will overwrite pre-existing results. 
            If False, will check whether a successfull calculation exists. If
            one does, then does not execute. 
        :param kwargs:
            Any attribute of the VASP instance can be overidden for
            the duration of this call by passing it as keyword argument.  

        :return: Yields an extractor object if a prior successful run exists.
                 Otherwise, yields a tuple object for executing an external
                 program.

        :note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        :raise RuntimeError: when computations do not complete.
        :raise IOError: when outdir exists but is not a directory.
    """ 
    from os.path import abspath
    from os import getcwd
    from ..process.program import ProgramProcess
    from ..misc import Changedir
    from .. import crystal_program

    # check for pre-existing and successfull run.
    if not overwrite:
      extract = self.Extract(outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return
    
    if outdir == None: outdir = getcwd()
    if workdir == None: workdir = self.guess_workdir(outdir)

    outdir = abspath(outdir)
    workdir = abspath(workdir)
    with Changedir(workdir) as tmpdir: 

      # writes/copies files before launching.
      self.bringup(structure, workdir, restart=self.restart)

      # figure out the program to launch.
      program = self.program if self.program is not None else crystal_program
      if hasattr(program, '__call__'): program = program(self)

      # now creates the process, with a callback when finished.
      def onfinish(process, error): self.bringdown(structure, workdir, outdir)
      yield ProgramProcess( program, cmdline=['<', 'crystal.input'],
                            outdir=outdir, onfinish=onfinish,
                            stdout='crystal.out', stderr='crystal.err',
                            dompi=True )
    # yields final extraction object.
    yield Extract(outdir)

  def __call__( self, structure, outdir=None, workdir=None, comm=None,
                overwrite=False, **kwargs):
    """ Calls CRYSTAL program. """
    result = None
    for program in self.iter( structure, outdir=outdir, workdir=workdir,
                              comm=comm, overwrite=overwrite, **kwargs ):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False):
        result = program
        continue
      # otherwise, it should yield a Program tuple to execute.
      program.start(comm)
      program.wait()
    # Last yield should be an extraction object.
    if not result.success:
      raise RuntimeError("CRYSTAL failed to execute correctly.")
    return result

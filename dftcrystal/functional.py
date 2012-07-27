__docformat__ = "restructuredtext en"
from ..functools import stateless, assign_attributes
from .extract import Extract

class Functional(object):
  """ Wrapper for the CRYSTAL program. 
  
      This object provides a pythonic interface to the CRYSTAL_ program. It is
      modelled loosely after CRYSTAL_'s input:

        - The OPTGEOM keyword of the first code block can be accessed through
          the :py:class:`optgeom <lada.dftcrystal.optgeom.OptGeom>`
          :py:attr:`attribute <optgeom>`:
          
          .. code-block:: python

            # enable geometry optimization
            functional.optgeom.enabled = True
            # change some of the keywords
            functional.optgeom.maxcycle = True

          Geometry optimization must be explicitly enabled, as is done in the
          first line above. 
        - The interface to the second block of input (basis-functions) can be
          accessed through the :py:class:`basis
          <lada.dftcrystal.basis.BasisSet>` :py:attr:`attribute <basis>`. It
          allows to set the basis set itself, as well as keywords specific to
          the basis set:

          .. code-block:: python
            
            functional.basis['H'] = [Shell('s', a0=(15.0, 1.0), a1=(8.0, 1.0)), ...]
            functional.basis.ghosts = [1, 4, 6], False

        - The third input block (Hamiltonian and miscellaneous) can be accessed
          directly through the functional, or, alternatively, *via* the
          :py:class:`scf <lada.dftcrystal.electronic.Electronic>`
          :py:attr:`attribute <scf>`.
        
          .. code-block:: python
        
            functional.scf.dft.b3lyp = True
            functional.dft.b3lyp = True
        
            functional.tolinteg = [8] * 4 + [14]
        
          The first two lines are exactly equivalent. The third line could also
          be written as ``functional.scf.tolinteg = ...``.

  """
  Extract = Extract
  """ Extraction class. """
  __ui_name__ = 'functional'
  """ Name used in user-friendly representation """
  def __init__(self, copy=None, program=None, **kwargs):
    """ Creates the crystal wrapper. """
    from .basis import BasisSet
    from .optgeom import OptGeom
    from .electronic import Electronic

    super(Functional, self).__init__()

    self.scf = Electronic()
    """ Holds scf/electronic keywords -- block 3. """
    self.basis   = BasisSet()
    """ Holds definition of basis functions -- block 2. """
    self.optgeom = OptGeom()
    """ Holds definition of geometry optimization -- part of block 1. """
    self.title   = None
    """ Title of the calculation. 
    
        Overriden by the name of the input structure, if it exists.
    """
    self.program = program
    """ Path to crystal program.

        If this attribute is None, then :py:data:`~lada.crystal_program` is
        used.
    """ 
    self.restart = None
    """ Place holder. """

  def __getattr__(self, name):
    """ Pushes scf stuff into instance namespace. """
    from ..error import ValueError
    if name in self.scf._crysinput: return getattr(self.scf, name)
    raise ValueError('Unknown attribute {0}.'.format(name))
  def __setattr__(self, name, value):
    """ Pushes scf stuff into instance namespace. """
    from .input import Keyword
    if isinstance(value, Keyword): 
      if name in ['scf', 'basis', 'optgeom']:
        return super(Functional, self).__setattr__(name, value)
      setattr(self.scf, name, value)
    elif name in self.scf._crysinput:
      setattr(self.scf, name, value)
    else: super(Functional, self).__setattr__(name, value)
  def __delattr__(self, name):
    """ Deletes attributes. """
  def __dir__(self):
    """ List of attributes and members """
    return list( set(self.__dict__.iterkeys()) | set(dir(self.__class__))      \
                 | set(self.scf._crysinput.iterkeys()) )

  def add_keyword(self, name, value=None):
    """ Passes on to :py:attr:`~Functional.scf` """
    return self.scf.add_keyword(name, value)

  def read_input(self, tree, owner=None):
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
      self.optgeom.read_input(tree[starter]['OPTGEOM'], owner=self)

    # read basis set
    if 'BASISSET' in tree.keys(): 
      self.basis.read_input(tree['BASISSET'], owner=self)

    # read hamiltonian stuff.
    if 'HAMSCF' in tree.keys():  
      self.scf.read_input(tree['HAMSCF'], owner=self)

  def print_input(self, **kwargs):
    """ Dumps CRYSTAL input to string. """

    # create optgeom part first, since it needs be inserted in the structure
    # bit. Split along lines and remove empty lines at end.
    # if empty, then make it empty.
    optgeom = self.optgeom.print_input(**kwargs)
    if optgeom is not None:
      optgeom = optgeom.rstrip().split('\n')
      while len(optgeom[-1].rstrip().lstrip()) == 0: optgeom.pop(-1)
      if len(optgeom) == 2: optgeom = []
    else: optgeom = []

    result = ''
    if 'structure' in kwargs:
      structure = kwargs['structure']
      # insert name of structure as title.
      if hasattr(structure, 'name'):
        result += str(structure.name).rstrip().lstrip() + '\n'
      elif getattr(self, 'title', None) is not None:
        result += str(self.title).rstrip().lstrip()
      result = result.rstrip()
      if len(result) == 0 or result[-1] != '\n': result += '\n'
       
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
    if len(result): result += '\n'
    result += self.basis.print_input(**kwargs)

    # add scf block
    result = result.rstrip()
    if len(result): result += '\n'
    result += self.scf.print_input(**kwargs)

    # end input and return
    result = result.rstrip()
    if len(result): result += '\n'
    return result + 'END\n'

  def guess_workdir(self, outdir):
    """ Tries and guess working directory. """
    from os import environ, getpid
    from os.path import join
    from datetime import datetime
    return join( environ.get('PBS_TMPDIR', outdir),
                 '{0!s}.{1}'.format(datetime.today(), getpid())\
                            .replace(' ', '_') )

  def bringup(self, structure, workdir, restart):
    """ Creates file environment for run. """
    from os.path import join
    from ..misc import copyfile, Changedir
    from ..error import ValueError
    from .. import CRYSTAL_filenames

    # sanity check
    if len(self.basis) == 0:
      raise ValueError('Basis is empty, cannot run CRYSTAL.')
      
    with Changedir(workdir) as cwd:
      # first copies file from current working directory
      if restart is not None: 
        for key, value in CRYSTAL_filenames.iteritems():
          copyfile( value.format('crystal'), key, nocopyempty=True,
                    symlink=False, nothrow="never" )
      # then copy files from restart.
      if restart is not None:
        for key, value in CRYSTAL_filenames.iteritems():
          copyfile( join(restart.directory, value.format('crystal')), 
                    key, nocopyempty=True, symlink=False, 
                    nothrow="never" )

      # then creates input file.
      string = self.print_input(crystal=self, structure=structure)
      with open('crystal.d12', 'w') as file: file.write(string)

  def bringdown(self, structure, workdir, outdir):
    """ Copies files back to output directory. 
    
        Cats input intO output.
	Removes workdir if different from outdir **and** run was successfull.
    """
    from os.path import join, samefile
    from shutil import rmtree
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
        file.write(self.__repr__(defaults=False))
        file.write('\n{0} END {1} {0}\n'.format(header, 'FUNCTIONAL'))
      if len([0 for filename in iglob(join(workdir, 'ERROR.*'))]):
        with open('crystal.err', 'w') as out:
          for filename in iglob(join(workdir, 'ERROR.*')):
            with open(filename, 'r') as file: out.write(file.read() + '\n')
    
    if Extract(outdir).success and not samefile(outdir, workdir):
      rmtree(workdir)
      

  
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
      dompi = comm is not None
      if dompi:
        from ..misc import copyfile
        copyfile('crystal.d12', 'INPUT')

      # figure out the program to launch.
      program = self.program if self.program is not None else crystal_program
      if hasattr(program, '__call__'):
        from inspect import getargspec
        args = getargspec(program)
        if 'comm' not in args.args and args.kwargs is None:
              program = program(self)
        else: program = program(self, comm=comm)

      # now creates the process, with a callback when finished.
      def onfinish(process, error): self.bringdown(structure, workdir, outdir)
      yield ProgramProcess( program, outdir=workdir, onfinish=onfinish,
                            stdout=None if dompi else 'crystal.out', 
                            stderr='crystal.out' if dompi else 'crystal.err',
                            stdin=None if dompi else 'crystal.d12', 
                            dompi=dompi )
    # yields final extraction object.
    yield Extract(outdir)

  def __call__( self, structure, outdir=None, workdir=None, comm=None,         \
                overwrite=False, **kwargs):
    for program in self.iter( structure, outdir=outdir, workdir=workdir,
                              comm=comm, overwrite=overwrite, **kwargs ):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False): continue
      # otherwise, it should yield a Program tuple to execute.
      program.start(comm)
      program.wait()
    # Last yield should be an extraction object.
    if not program.success:
      raise RuntimeError("CRYSTAL failed to execute correctly.")
    return program
  __call__.__doc__ = iter.__doc__
 
  def __repr__(self, defaults=True, name=None):
    """ Returns representation of this instance """
    from ..functools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    from ..functools.uirepr import template_ui_repr

    results = template_ui_repr(self, imports, name, defaults, ['scf'])
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    scf = self.scf.__ui_repr__(imports, name, getattr(defaults, 'scf', None))
    results.update(scf)
    return results

  def __deepcopy__(self, memo):
    from copy import deepcopy
    result = self.__class__()
    result.__dict__ = deepcopy(self.__dict__)
    return result

  def __getstate__(self): return self.__dict__
  def __setstate__(self, value):
    self.__dict__.update(value.copy())

  def copy(self): 
    from copy import deepcopy
    return deepcopy(self)

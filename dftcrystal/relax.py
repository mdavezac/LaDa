__docformat__ = "restructuredtext en"
__all__ = ['relax', 'iter_relax', 'Relax', 'RelaxExtract']
from ..tools.makeclass import makeclass, makefunc
from .functional import Functional
from .extract import Extract, MassExtract

class RelaxExtract(Extract):
  """ Extractor class for vasp relaxations. """
  class IntermediateMassExtract(MassExtract):
    """ Focuses on intermediate steps. """
    def __iter_alljobs__(self):
      """ Goes through all directories with an OUTVAR. """
      from glob import iglob
      from os.path import relpath, join, exists

      for dir in iglob(join(join(self.rootpath, 'relax'), '*/')):
        if not exists(join(self.rootpath, join(dir, 'crystal.out'))): continue
        try: result = Extract(dir[:-1])
        except: continue
        yield join('/', relpath(dir[:-1], self.rootpath)), result
        
  @property
  def details(self):
    """ Intermediate steps. """
    if '_details' not in self.__dict__:
      self.__dict__['_details'] = self.IntermediateMassExtract(self.directory)
    return self._details
  
  def iterfiles(self, **kwargs):
    """ Iterates over input/output files. """
    from itertools import chain
    for file in chain( super(RelaxExtract, self).iterfiles(**kwargs),
                       self.details.iterfiles(**kwargs) ): yield file

def iter_relax(self, structure=None, outdir=None, maxiter=30, **kwargs):
  """ Performs relaxations until convergence is reached.
  
      A final static calculation is performed at the end of the run.
  """
  from os.path import join
  from os import getcwd
  from ..error import ValueError, ExternalRunFailed

  self = self.copy()
  if maxiter <= 0:
    raise ValueError('Maximum number of iteration cannot be less than 0.')
  if outdir is None: outdir = getcwd()
  self.optgeom.enabled = True
  if self.optgeom.maxcycle is None: 
    self.optgeom.maxcycle = 10
  elif self.optgeom.maxcycle < 0: 
    self.optgeom.maxcycle = 10
   
  
  iteration, restart = 0, kwargs.pop('restart', None)

  while iteration < maxiter:
    # performs calculation
    outpath = join(join(outdir, 'relax'), str(iteration))
    for extract in self.iter( structure, outdir=outpath, restart=restart,
                              **kwargs ):
      yield extract

    # check for success
    if not extract.success: 
      raise ExternalRunFailed( 'CRYSTAL run did not succeed in {0}.'           \
                               .format(outpath) )

    # otherwise, go to next iteration
    structure = extract.crystal
    restart   = extract
    iteration += 1

    # check for convergence
    if extract.optgeom_convergence is True                                     \
       and extract.optgeom_iterations < self.optgeom.maxcycle: break

  # perform static calculation
  self.optgeom.enabled = False
  for extract in self.iter( structure, outdir=outdir, restart=restart, 
                            **kwargs ):
    yield extract

  yield iter_relax.Extract(outdir)     

iter_relax.Extract = RelaxExtract
""" Extraction object for relaxation meta-functional. """

Relax = makeclass( 'Relax', Functional, iter_relax, None,
                   module='lada.dftcrystal.relax',
		   doc='Functional form of the '                               \
                       ':py:class:`lada.dftcrystal.relax.iter_relax` method.' )
relax = makefunc('relax', iter_relax, module='lada.dftcrystal.relax')

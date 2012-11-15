from ..tools import stateless, assign_attributes
from .vff import Vff
from .extract import Extract as ExtractVFF

class Functional(Vff): 
  Extract = ExtractVFF
  """ Extraction object for Vff. """
  def __init__( self, relax=True, method='BFGS', tol=1e-8, maxiter=100,
                verbose=True, copy=None, options=None ): 
    super(Functional, self).__init__()

    self._parameters = {}
    """ Holds vff parameters. """
    self.relax = relax
    """ Whether to relax the structure """
    self.method = method
    """ Type of method used to relax the structure. 
    
        .. see:: 
         
          `scipy.optimize.minimize`__'s method argument.
          
        .. __: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    """
    self.tol = tol
    """ Convergence criteria. """
    self.maxiter = maxiter
    """ Maximum number of iterations. """
    self.verbose = verbose
    """ Whether minimization should be verbose. """
    self.options = options
    """ Additional options for the chosen minimizer. """

    if copy is not None: self.__dict__.update(copy.__dict__)


  def _is_static(self, **kwargs):
    """ True if calculation is static. """
    relax = kwargs.get('relax', self.relax)
    if relax is False or relax is None: return True
    if not isinstance(relax, str): return False
    return relax.lower() == 'static'

  @stateless
  @assign_attributes(ignore=['overwrite', 'comm'])
  def __call__(self, structure, outdir=None, overwrite=False, **kwargs):
    """ Evaluates energy and forces on a structure. """
    from datetime import datetime
    from os.path import join
    from scipy.optimize import minimize
    from .extract import Extract as ExtractVFF
    from ..misc import Redirect

    if ExtractVFF(outdir).success and not overwrite: return ExtractVFF(outdir)

    header = ''.join(['#']*20)
    with open(join(outdir, 'vff.out'), 'w') as file:
      file.write('Start date: {0!s}\n'.format(datetime.today()))
      file.write('{0} {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
      file.write( 'from {0.__class__.__module__} '                             \
                  'import {0.__class__.__name__}\n'.format(structure) )
      string = repr(structure).replace('\n', '\n            ')
      file.write('structure = ' + string + '\n')
      file.write('{0} END {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
      file.write('{0} {1} {0}\n'.format(header, 'FUNCTIONAL'))
      file.write(self.__repr__(defaults=False) + '\n')
      file.write('{0} END {1} {0}\n\n'.format(header, 'FUNCTIONAL'))
    minimization, result = None, None
    try: 
      if self._is_static(**kwargs):
        result = super(Functional, self).__call__(structure)
      else: 
        funcs = self._getfuncs_relaxall(structure)
        options = {} if self.options is None else self.options.copy()
        options['disp'] = self.verbose
        options['maxiter'] = self.maxiter
        with Redirect(join(outdir, 'vff.out'), ['out', 'err'], True) as f:
          minimization = minimize( funcs.energy, jac=funcs.jacobian, 
                                   x0=funcs.x0, tol=self.tol,
                                   method=self.method, options=options )
        result = funcs.structure(minimization.x)
        # this will reconstruct the tree and might fail. 
        result = super(Functional, self).__call__(result)
    finally: 
      with open(join(outdir, 'vff.out'), 'a') as file:
        if minimization is not None:
          file.write('{0} {1} {0}\n'.format(header, 'MINIMIZATION'))
          file.write(repr(minimization) + '\n')
          file.write('{0} END {1} {0}\n\n'.format(header, 'MINIMIZATION'))
        if result is not None:
          file.write('{0} {1} {0}\n'.format(header, 'STRUCTURE'))
          file.write( 'from {0.__class__.__module__} '                         \
                      'import {0.__class__.__name__}\n'.format(result) )
          string = repr(result).replace('\n', '\n            ')
          file.write('structure = ' + string + '\n')
          file.write('{0} END {1} {0}\n'.format(header, 'STRUCTURE'))
        file.write('End date: {0!s}\n'.format(datetime.today()))
    return ExtractVFF(outdir)


  def _getfuncs_relaxall(self, structure):
    """ Functions when relaxing all degrees of freedom. """
    from collections import namedtuple
    from numpy import dot, array, zeros
    from numpy.linalg import inv, det
    from quantities import angstrom
    from . import build_tree

    structure = structure.copy()
    cell0, pos0 = structure.cell.copy(), structure[0].pos
    tree = build_tree(structure)
    scale = float(structure.scale.rescale(angstrom))

    def xtostrain(x0):
      return array([[x0[0] + 1e0, x0[1], x0[2]],
                    [x0[1], 1e0 + x0[3], x0[4]],
                    [x0[2], x0[4], x0[5]+1e0]])
      
    def update_structure(x0, strain):
      structure.cell = dot(strain, cell0)
      structure[0].pos = dot(strain, pos0)
      for i, atom in enumerate(structure[1:]):
        atom.pos = dot(structure.cell, x0[i*3:3+i*3])

    def make_structure(x0):
      strain = xtostrain(x0)
      update_structure(x0[6:], strain)
      return structure

    def energy(x0):
      strain = xtostrain(x0)
      update_structure(x0[6:], strain)
      return self.energy(structure, _tree=tree).magnitude

    def jacobian(x0):
      strain = xtostrain(x0)
      update_structure(x0[6:], strain)

      stress, forces = self.jacobian(structure, _tree=tree)
      stress *= -det(structure.scale * structure.cell)

      stress = dot(stress, inv(strain))
      result = zeros(3+len(structure)*3, dtype='float64')
      result[0] = stress[0,0]
      result[1] = 2e0*stress[0,1]
      result[2] = 2e0*stress[0, 2]
      result[3] = stress[1,1]
      result[4] = 2e0*stress[1,2]
      result[5] = stress[2, 2]
      result[6:] = dot(structure.cell.T*scale, forces[1:].T).T.flatten()
      return result

    x = zeros(3+len(structure)*3, dtype='float64')
    frac = inv(structure.cell)
    for i, atom in enumerate(structure[1:]): x[6+3*i:9+3*i] = dot(frac, atom.pos)

    Functions = namedtuple('Functions', ['x0', 'jacobian', 'energy', 'structure'])
    return Functions(x, jacobian, energy, make_structure)

del stateless
del assign_attributes
del Vff
del ExtractVFF

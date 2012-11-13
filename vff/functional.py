from .vff import Vff
class Functional(Vff): 
  def __init__(self, relax=True, method='BFGS', tol=1e-8, maxiter=50): 
    super(Functional, self).__init__()

    self._parameters = {}
    """ Holds vff parameters. """
    self.relax = relax
    """ Whether to relax the structure """
    self.methods = method
    """ Type of method used to relax the structure. 
    
        .. see:: 
         
          `scipy.optimize.minimize`__'s method argument.
          
        .. __: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    """
    self.tol = tol
    """ Convergence criteria. """
    self.maxiter = maxiter
    """ Maximum number of iterations. """


  def _is_static(self, **kwargs):
    """ True if calculation is static. """
    relax = kwargs.get('relax', self.relax)
    if relax is False or relax is None: return True
    if not isinstance(relax, str): return False
    return relax.lower() == 'static'


  def __call__(self, structure, **kwargs):
    """ Evaluates energy and forces on a structure. """
    if self._is_static(**kwargs):
      result = super(self, Functional).__init__(structure)

    else: 
      result = self._relax_all(structure)

    return result


  def _relax_all(self, structure):
    from numpy import dot, array, zeros
    from numpy.linalg import inv, det
    from scipy.optimize import fmin_bfgs as minimize
    from . import build_tree

    structure = structure.copy()
    cell0 = structure.cell.copy()
    tree = build_tree(structure)

    def xtostrain(x0):
      return array([[x0[0] + 1e0, x0[1], x0[2]],
                    [x0[1], 1e0 + x0[3], x0[4]],
                    [x0[2], x0[4], x0[5]+1e0]])
      
    def update_structure(x0, strain):
      structure.cell = dot(strain, cell0)
      for i, atom in enumerate(structure):
        atom.pos = dot(structure.cell, x0[i*3:3+i*3])

    def gradient_tox(stress, forces, strain):
      from numpy import dot
      from numpy.linalg import inv

      stress = dot(stress, inv(strain))
      result = stress[0].tolist() + stress[1,1:].tolist() + [stress[2,2]]
      result += dot(inv(structure.cell), forces.T).T.flatten().tolist()
      return array(result)


    def energy(x0):
      strain = xtostrain(x0)
      update_structure(x0[6:], strain)
      return self.energy(structure, _tree=tree).magnitude

    def jacobian(x0):
      strain = xtostrain(x0)
      update_structure(x0[6:], strain)

      stress, forces = self.jacobian(structure, _tree=tree)
      stress *= -det(structure.scale * structure.cell)
      return gradient_tox(stress, forces, strain)

    x = zeros(6+len(structure)*3, dtype='float64')
    x[:6] = 0e0
    frac = inv(structure.cell)
    for i, atom in enumerate(structure): x[6+3*i:9+3*i] = dot(frac, atom.pos)

    return minimize( energy, fprime=jacobian, x0=x, 
                     gtol=self.tol, maxiter=self.maxiter )

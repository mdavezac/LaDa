__docformat__ = "restructuredtext en"
__all__ = ['Electronic']
from input import AttrBlock, TypedKeyword, BoolKeyword

class Electronic(AttrBlock):
  """ DFT attribute block. """ 
  def __init__(self):
    """ Creates the scf attribute block. """
    from .hamiltonian import Dft
    super(Electronic, self).__init__()
    self.maxcycle = TypedKeyword(type=int)
    """ Maximum number of electronic minimization steps """
    self.tolinteg = TypedKeyword(type=[int]*5)
    """ Integration truncation criteria """
    self.toldep   = TypedKeyword(type=int)
    """ Density matrix convergence criteria """
    self.tolpseud = TypedKeyword(type=int)
    """ Pseudopotential truncation criteria """
    self.toldee   = TypedKeyword(type=int)
    """ Total energy convergence criteria """
    self.testpdim = BoolKeyword()
    """ Stop after processing input and performin symmetry analysis """
    self.test     = BoolKeyword()
    """ Stop after printing ressource requirement """
    self.symadapt = BoolKeyword()
    """ Symmetry adapted bloch wavefunctions """
    self.savewf   = BoolKeyword()
    """ Save wavefunctions to disk """
    self.dft      = Dft()
    """ Holds definition of the DFT functional itself """

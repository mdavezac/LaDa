from parameters import Incar


class Vasp(object):

  program = "vasp.mpi"
  def __init__(self):

    # incar parameters.
    self.incar = Incar()

  def _prerun(self):

    for param in self.incar:
      print param

  def __call__(self, structure):

    self.sytem = structure
    self._prerun()


if __name__ == "__main__":
  from lada import crystal, atat
  vasp = Vasp() 
  vasp.species = [("Al", "~/AtomicPotentials/pseudos/K_s")]
  structure = crystal.sStructure()
  structure.scale = 1e0
  structure.cell = atat.rMatrix3d([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
  structure.atoms.append( crystal.StrAtom() )
  structure.atoms[0].pos = atat.rVector3d(0,0,0)
  structure.atoms[0].type = "Al"

  print structure

  vasp(structure)








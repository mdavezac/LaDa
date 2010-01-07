from incar import Incar
from kpoints import Density


class Vasp(Incar):

  program = "vasp.mpi"

  def __init__(self):
    Incar.__init__(self) 

    self.indir = ""
    self.kpoints = Density()

  def _prerun(self):

    for param in self:
      print param.incar_string(self)

    print self.kpoints(self)

  def __call__(self, structure):

    self.system = structure
    self._prerun()


if __name__ == "__main__":
  from lada import crystal, atat
  from specie import Specie
  
  vasp = Vasp() 
  vasp.species = [Specie("Al", "~/AtomicPotentials/pseudos/K_s")]
  vasp.fftgrid.value = (10,10,10)
  structure = crystal.sStructure()
  structure.scale = 1e0
  structure.cell = atat.rMatrix3d([[2,0.5,0.5],[0.0,0,0.5],[0.0,0.5,0]])
  structure.atoms.append( crystal.StrAtom() )
  structure.atoms[0].pos = atat.rVector3d(0,0,0)
  structure.atoms[0].type = "Al"

  print structure

  vasp(structure)








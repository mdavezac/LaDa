""" Module to extract esca and vff ouput. """

from ..vff import Extract as VffExtract, _get_script_text
from ..opt.decorators import broadcast_result, make_cached

class Extract(VffExtract):
  def __init__(self, directory = None, comm = None, escan = None):

    if escan == None: escan = Escan()
    super(Extract, self).__init__(directory, comm = comm, vff = escan.vff)
  
  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for VFF success.
        
        At this point, checks for files and 
    """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    if not exists(path): return False

    with open(path, "r") as file:
      for line in file:
        if line.find("# Computed in:") != -1: return True
    return False

  @property
  @make_cached
  def escan(self):
    """ Gets escan functional from self.L{OUTCAR}. """
    from os.path import exists, join
    from .. import Escan, localH, nonlocalH, soH, AtomicPotential
    
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (file))

    @broadcast_result(attr=True, which=0)
    def get_functional(this):
      with open(path, "r") as file: return _get_script_text(file, "Escan")
    local_dict = { "lattice": self.lattice, "minimizer": self.minimizer,\
                   "vff": self.vff, "Escan": Escan, "localH": localH,\
                   "nonlocalH": nonlocalH, "soH": soH, "AtomicPotential":AtomicPotential}
    exec get_vff(self) in globals(), local_dict
    return local_dict["escan"]

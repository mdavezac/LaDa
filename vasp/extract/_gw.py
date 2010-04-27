""" Extracts VASP-GW output """
from . import _ExtractImpl
from .decorators import bound_broadcast_result 
from ...opt.decorators import make_cached
class _ExtractGWImpl(_ExtractImpl):
  """ Implementation of GW extractor. """
  def __init__(self, *args, **kwargs):
    """ Same as L{_ExtractImpl} """
    super(_ExtractGWImpl, self).__init__(*args, **kwargs)

  def _get_eigocc(self,which):
    """ Implementation of _get_eigenvalues and _get_occupations """
    import re 
    from os.path import exists, join
    from numpy import array

    path = self.OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      found = re.compile(r"""k-point\s+(\d+)\s*:\s*(\S+)\s+(\S+)\s+(\S+)$""")
      in_kpoint = -1
      kp_result = []
      for line in file:
        if in_kpoint > 1: 
          data = line.split()
          if len(data) == 8: kp_result.append(float(data[which]))
          else: 
            result.append(kp_result)
            kp_result = []
            in_kpoint = -1
        elif in_kpoint >= 0 and in_kpoint < 2: in_kpoint += 1
        else:
          match = found.search(line)
          if match != None:  
            if int(match.group(1)) == 1: result = []
            in_kpoint = 0
    return array(result, dtype="float64")

  @make_cached
  @bound_broadcast_result
  def _get_dft_eigenvalues(self):
    """ Returns DFT eigenvalues 

        @return: a two-dimension numpy nxm array of eigenvalues, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(1)
  
  @make_cached
  @bound_broadcast_result
  def _get_qp_eigenvalues(self):
    """ Returns Quasi-Particle eigenvalues 

        @return: a two-dimension numpy nxm array of eigenvalues, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(2)
  
  @make_cached
  @bound_broadcast_result
  def _get_self_energies(self):
    """ Returns self-energies for each eigenvalue

        @return: a two-dimension numpy nxm array of eigenvalues, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(3)
  
  @make_cached
  @bound_broadcast_result
  def _get_occupations(self):
    """ Returns occupation for each eigenvalue

        @return: a two-dimension numpy nxm array of eigenvalues, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(7)


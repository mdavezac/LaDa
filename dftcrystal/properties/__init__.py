""" Subpackage grouping single-electron properties. """
__docformat__ = "restructuredtext en"
__all__ = ['Properties']
from ..input import AttrBlock

class Properties(AttrBlock):
  """ Wrap single-electron property calculations. """
  def __init__(self, input=None):
    from .keywords import NewK, Band, Rdfmwf
    from ..input import BoolKeyword
    super(Properties, self).__init__()
    self.rdfmwf = Rdfmwf()
    """ Whether to create binary wavefunction file. 
    
        - If None, LaDa figures out which of  'crystal.f98' and 'crystal.f9'
          file is the latest in input and output directories. Then figures out
          whether to include an RDFMWF tag or not.
        - If True, LaDa looks only for a 'crystal.f98' in the input and/or
          output directories, using the latest. If no file exists, then it is
          an error.
        - If False, LaDa checks that 'crystal.f9' exists in the input and/or
          output directories, using the latest. If no file exists, then it is
          an error.

        If no error occurs, the relevant files are copied to the working
        directory.
    """
    self.newk = NewK()
    """ Performs diagonalization on new k-point mesh. 

        Starting from the hamiltonian defined in input, performs a
        diagonalization on a new k-point mesh.
        The input is similar to :py:attr:`~lada.functional.shrink`.
        Additionaly, it is possible to specify printing options:

        >>> properties.newk.printing[66] = -5

        Following the same input scheme as
        :py:attr:`~lada.functional.setprint`.

        It is possible to set the k-point mesh directly, much as for
        :py:attr:`~lada.functional.shrink`:

        >>> properties.newk = 5, None

        However, both the ``recompute_fermi`` and the printing options
        ``printing`` must be set by direct access:

        >>> properties.newk.printing[5] = 0
        >>> properties.recompute_fermi = True
    """
    self.nosymada = BoolKeyword()
    """ No symmetry adapted bloch function. """
    self.band = Band()
    """ Performs calculation on a path of k-points.
    
        This attribute can be parameterized in a variety of way.

        - ``properties.band.title``: accepts a string which will be the title
          of the band-structure. This is mostly to reproduce original CRYSTAL_
          input.
        - ``properties.band.nsub``: accepts an integer which defines the total
          number of points along the path.
        - ``properties.band.iss``: see CRYSTAL_'s user guide. Not really
          relevant in python.
        - ``properties.band.minband``: Index of the minimum band for which to
          output results. Best leave this to zero since bands starting from
          zero are computed anyways. Hold-out from CRYSTAL_.
        - ``properties.band.maxband``: Index of the maximum band for which to
          ouput results.

        The path is held the ``lines`` attributes. It can be given as:

        >>> properties.band = [startA, endA], [startB, endB], ...
    """
    self.input = input

  def output_map(self, **kwargs):
    """ Returns map of crystal input. """
    from ...tools.input import Tree
    from ...misc import RelativePath
    if 'properties' not in kwargs: kwargs['properties'] = self
    if 'input' not in kwargs: kwargs['input'] = self.input
    if 'outdir' not in kwargs: kwargs['outdir'] = None
    kwargs['outdir'] = RelativePath(kwargs['outdir']).path 
    if 'workdir' not in kwargs: kwargs['workdir'] = None
    kwargs['workdir'] = RelativePath(kwargs['workdir']).path 
    root = Tree()
    # First add RDFWF
    AttrBlock._output_map(root, 'rdfmwf', self._input['rdfmwf'], **kwargs)
    # Then go on to symmetry adapted thingie
    AttrBlock._output_map(root, 'nosymada', self._input['nosymada'], **kwargs)
    # Move on to newk
    AttrBlock._output_map(root, 'newk', self._input['newk'], **kwargs)
    # Do other preliminaries.
    for key in ['pato', 'pban', 'pgeomw', 'pdide', 'pscf']:
      if key in self._input: 
        AttrBlock._output_map(root, key, self._input[key], **kwargs)
      elif key.upper() in self._input:
        AttrBlock._output_map(root, key.upper(), self._input[key], **kwargs)
    # Now do all others.
    prelims = set([ 'rdfmwf', 'nosymada', 'newk', 'pato', 
                    'pban', 'pgeomw', 'pdide', 'pscf' ]) 
    for key, value in self._input.iteritems():
      if key not in prelims and key.upper() not in prelims: 
        AttrBlock._output_map(root, key, value, **kwargs)
    return root
  
  def print_input(self, **kwargs):
    """ Prints input to string. """
    from ..input import print_input
    map = self.output_map(**kwargs)
    if len(map) == 0: return ""
    # Otherwise, everything is standard.
    return print_input(map).rstrip() + '\nEND\n'

  @property
  def input(self):
    """ Input calculation from which to start. 

        - If None, then an extraction object is created in the current
          directory.
        - If a string, then it is assumed to the path to the directory
          containing the input calculations.
        - If an extraction object, then that extraction object should point to
          the self-consistent calculations.
    """
    return self._input_calc
  @input.setter
  def input(self, value):
    from ..extract import Extract as CrystalExtract
    if value is None: self._input_calc = CrystalExtract()
    elif isinstance(value, str): self._input_calc = CrystalExtract(value)
    else: self._input_calc = value

  def nbelectrons(self, structure, input=None):
    if input is None: input = self.input
    species = [u.type for u in structure.eval()]
    result = 0
    for specie in set(species):
      result += sum(u.charge for u in input.functional.basis[specie])          \
                * species.count(specie)
    return result


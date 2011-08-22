""" Subpackage for gw calculations. """
__docformat__ = "restructuredtext en"
__all__ = ['GWFunctional']
from .functional import Functional as VaspFunctional

class GWFunctional(VaspFunctional):
  """ Performs GW calculations. """
  def __init__(self, vasp=None, empties=None, gwiterations=5, gwparams=None, **kwargs):
    """ Initializes GW functional.
    
    
    """
    from copy import deepcopy
    super(GWFunctional, self).__init__(**kwargs)
    self.empties = {} if empties == None else deepcopy(empties)
    """ Parameters for empty bands DFT calculations. """
    self.gwparams = { "algo": "gw", "lmaxfockae": 4, "nomega": 64, "precfock": False, 
                      "encutgw": 150, "encutlf": 150, "lrpa": True, "nelm": 1, "loptics": True,
                      "lpead": True }
    """ Parameters for actual GW calculations. """
    if gwparams != None: self.gwparams.update(gwparams)
    self.gwiterations = gwiterations
    """ Number of gwiterations to perform. """
    

  def generator(self, structure, outdir=None, comm=None, **kwargs):
    """ Performs all vasp calculation leading to and including GW calculations. """
    from os import getcwd
    from os.path import join
    from copy import deepcopy
    from ..opt import RelativeDirectory
    
    outdir = getcwd() if outdir == None else RelativeDirectory(outdir).path
    empties = deepcopy(self.empties)
    if 'empties' in kwargs: empties.update(kwargs.pop('empties'))
    gwparams = deepcopy(self.gwparams)
    if 'gwparams' in kwargs: gwparams.update(kwargs.pop('gwparams'))
    overwrite = kwargs.pop('overwrite', False)
    # first, creates functional with copied keywords.
    this = deepcopy(self)
    for key, value in kwargs.iteritems():
      assert hasattr(this, key), KeyError('{0} does not correspond to a vasp attribute.'.format(key))
      setattr(this, key, value)



    # Performs empty band calculation.
    # check for existence of empty bands calculations.
    empty_bands = self.Extract(join(outdir, join("GW", "empties")), comm=comm)
    if overwrite == False and not empty_bands.success:
      empties_func = VaspFunctional(vasp=this, **empties)
      # check that we are performing DFT calculation, otherwise gets default.
      if empties_func.algo not in ['gw', 'gw0', 'chi', 'scgw', 'scgw0']: 
        empties_func.algo = VaspFunctional().algo
      assert empties_func.relaxation == "static", ValueError("Cannot perform relaxation in GW.")
      # perform calculations.
      empty_bands = empties_func\
                    (
                      structure=structure,
                      outdir=join(outdir, join("GW", "empties")),
                      comm=comm, overwrite=overwrite, minversion=5
                    )
    yield empty_bands 
    assert empty_bands.success, RuntimeError("Empty band calculation was not successfull, aborting.")

    for iteration in xrange(self.gwiterations): 
      gwfunc = VaspFunctional(vasp=this, **gwparams)
      gwextract =  gwfunc\
                   (
                     structure=structure,
                     outdir=join(outdir, join("GW", str(iteration))),
                     comm=comm, overwrite=overwrite, minversion=5
                   )
      yield gwextract
      assert gwextract.success,\
             RuntimeError("GW iteration was not successfull, aborting.")

  def __call__(self, structure, outdir=None, comm=None, **kwargs):
    """ Performs all steps of GW calculations. """
    for e in self.generator(structure, outdir, comm, **kwargs): continue
    return e


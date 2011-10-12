""" Subpackage for gw calculations. """
__docformat__ = "restructuredtext en"
__all__ = ['GWFunctional']
from .functional import Functional as VaspFunctional

class GWFunctional(VaspFunctional):
  """ Performs GW calculations. """
  def __init__(self, gwiterations=5, empties=None, gwparams=None, vasp=None, **kwargs):
    """ Initializes GW functional.
    
        :Parameters:
          gwiterations : int 
            Number of separate GW optimization to perform.
          empties : dict
	    Parameters specific to computations of empty bands. By default,
	    contains only ``loptics = true``. This value is always present
	    unless it is specifically overriden. Finally, parameters in
            ``kwargs`` take precedence.
          gwparams : dict
	    Parameters specific to gw calculations. These values are superseded
            by those in ``kwargs``. 
	    By default, contains the following parameters:

            - algo: "gw"
            - lmaxfockae: 4
            - nomega: 64
            - precfock: 'fast'
            - encutgw: 150
            - encutlf: 150
            - lrpa: True
            - nelm: 1
            - loptics: True
            - lpead: True

	    These default value are always present unless specifically
	    overriden. Finally, parameters in ``kwargs`` take precedence.
      The "npar" parameter is set to the number of processes (comm.size) when
      performing GW calculations. This cannot be overriden at this time since
      VASP bails out on any other value.
          vasp : None or vasp functional
	    If ``vasp`` is present, it will be cloned into this instance, e.g.
	    the parameters of this instance will be the same as those of
	    ``vasp``. Extra keyword arguments (from ``kwargs``), from
	    ``empties``, and from ``gwparams`` take precedence over parameters
            from the cloned functional.
          kwargs 
	    keyword arguments to set other vasp parameters for all calculations
	    (empty bands and gw). These values supersede all others, including
            those from ``vasp``, ``empties``, and ``gwparams``.
    """
    from copy import deepcopy

    super(GWFunctional, self).__init__(**kwargs)

    self.empties = {'loptics': True}
    """ Parameters for empty bands DFT calculations. """
    if empties is not None: self.empties.update(empties)
    for key in deepcopy(self.empties): 
      if key in kwargs: del self.empties[key]
    self.gwparams = { "algo": "gw", "lmaxfockae": 4, "nomega": 64, "precfock": 'fast',
                      "encutgw": 150, "encutlf": 150, "lrpa": True, "nelm": 1, "loptics": True,
                      "lpead": True }
    if gwparams is not None: self.gwparams.update(gwparams)
    for key in deepcopy(self.gwparams): 
      if key in kwargs: del self.gwparams[key]
    """ Parameters for actual GW calculations. """
    self.gwiterations = gwiterations
    """ Number of gwiterations to perform. """
    

  def generator(self, structure, outdir=None, comm=None, **kwargs):
    """ Performs all vasp calculation leading to and including GW calculations. """
    from os import getcwd
    from os.path import join
    from copy import deepcopy
    from ..opt import RelativeDirectory
    
    outdir = getcwd() if outdir is None else RelativeDirectory(outdir).path
    empties = deepcopy(self.empties)
    if 'empties' in kwargs: empties.update(kwargs.pop('empties'))
    gwparams = deepcopy(self.gwparams)
    if 'gwparams' in kwargs: gwparams.update(kwargs.pop('gwparams'))
    overwrite = kwargs.pop('overwrite', False)
    if 'minversion' in kwargs:
      assert kwargs.pop('minversion') >= 5, \
             ValueError("Requested vasp version < 5 for GW calculations.")
    if "npar" in kwargs: 
      empties["npar"] = kwargs.pop("npar")
    gwparams["npar"] = comm.size

    # Performs empty band calculation.
    # check for existence of empty bands calculations.
    empty_bands = self.Extract(join(outdir, join("GW", "empties")), comm=comm)
    if overwrite == False and not empty_bands.success:
      empties_func = VaspFunctional(vasp=self, **empties)
      # check that we are performing DFT calculation, otherwise gets default.
      if empties_func.algo in ['gw', 'gw0', 'chi', 'scgw', 'scgw0']: 
        empties_func.algo = VaspFunctional().algo
      assert empties_func.relaxation == "static", ValueError("Cannot perform relaxation in GW.")
      assert empties_func.loptics, ValueError("Cannot perform GW calculations without loptics=True.")
      # perform calculations.
      empty_bands = empties_func\
                    (
                      structure=structure,
                      outdir=join(outdir, join("GW", "empties")),
                      comm=comm, overwrite=overwrite, minversion=5,
                      **kwargs
                    )
    yield empty_bands 
    assert empty_bands.success, RuntimeError("Empty band calculation was not successfull, aborting.")

    kwargs['restart'] = empty_bands
    for iteration in xrange(self.gwiterations): 
      gwfunc = VaspFunctional(vasp=self, **gwparams)
      assert gwfunc.npar == comm.size
      assert gwfunc.relaxation == "static", ValueError("Cannot perform relaxation in GW.")
      assert gwfunc.loptics, ValueError("Cannot perform GW calculations without loptics=True.")
      assert gwfunc.algo in ['gw', 'gw0', 'scgw', 'scgw0', 'chi'], \
             ValueError("Requested gw calculation without gw algorithm.")
      gwextract =  gwfunc\
                   (
                     structure=structure,
                     outdir=join(outdir, join("GW", str(iteration))),
                     comm=comm, overwrite=overwrite, minversion=5,
                     **kwargs
                   )
      yield gwextract
      assert gwextract.success,\
             RuntimeError("GW iteration was not successfull, aborting.")
      kwargs['restart'] = gwextract

  def __call__(self, structure, outdir=None, comm=None, **kwargs):
    """ Performs all steps of GW calculations. """
    for e in self.generator(structure, outdir, comm, **kwargs): continue
    return e


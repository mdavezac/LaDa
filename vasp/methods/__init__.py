#! /usr/bin/python
""" Module with different methods for vasp 

    Most of the generators in this module will not recompute the result if an
    appropriate output directory is found. In other words, this function can be
    called once to perform calculation, and then to perform output extractions.
    >>> # First compute results (if not yet computed).
    >>> for extract in relaxation(structure, vasp, "relaxation"): pass
    >>> # do something else
    >>> # Then  analyze relaxation steps.
    >>> for extract in relaxation(structure, "relaxation", vasp):
    >>>   print "Where: ", extract.indir
    >>>   print extract.total_energy
    The output extraction object will be the output the vasp callable.
"""
def relaxation( structure, vasp, outdir="relaxation", repat = [], tolerance = 1e-3, \
                relaxation="volume ionic cellshape", **kwargs ):
  """ Performs a vasp relaxation

      This is a generator which yields a output extraction object
      after each vasp calculation:
      >>> for extract in relaxation(structure, vasp, "outdir")
      >>>   print extract.total_energy
      Note that this function will not recompute the result if an appropriate
      output directory is found. 
      @note: This functor is stateless as long as self and structure can be
             deepcopied correctly.  

      @param structure: A structure to relax.
      @type structure: L{lada.crystal.sStructure} or L{lada.crystal.Structure}.
      @param vasp: The vasp functional.
      @type vasp: L{Vasp}
      @param outdir: Directory where to repatriate files. Default = "relaxation"
      @type outdir: str
      @param repat: File to repatriate, other than L{files.minimal}. Default: [].
      @type repat: list or set
      @param tolerance: Total energy convergence criteria. Default: 1e-3. 
      @type tolerance: float
      @param relaxation: Which degrees of freedom to relax. Default: \"volume
        ionic cellshape\". @see L{incar.Incar.relaxation}.
      @type relaxation: string. 
  """
  from copy import deepcopy
  from math import fabs as abs
  from os.path import join, exists
  from os import remove
  from .. import files

  # make this function stateless.
  vasp = deepcopy(vasp)
  structure = deepcopy(structure)
  repat = set(repat).union(files.minimal)
  tolerance = float(tolerance)
  vasp.relaxation = str(relaxation)

  # number of restarts.
  nb_steps = 0 

  # we may want to include restart file to speed up calculations. However,
  # these files should be deleted afterwards.
  other_repat = []
  for file in files.restart: 
    if file not in repat: other_repat.append(file)
  repat.union(other_repat)

  # performs initial calculation.
  outdirs = ["%s/step_%i" % (outdir, nb_steps)]
  output = vasp(structure, outdirs[-1], repat, **kwargs)
  # yields output for whatnot
  yield output

  structure = crystal.sStructure(output.structure)
  # makes sure we don't accidentally converge to early.
  oldenergy = -structure.energy

  tol = tolerance * float( len(structure.atoms) )
  while( abs(oldenergy - structure.energy) > tol ):
    # restart + new directory
    vasp.indir = outdirs[-1]
    nb_steps += 1

    # launch calculations 
    outdirs.append("%s/step_%i" % (outdir, nb_steps))
    output = vasp(structure, outdirs[-1], repat, **kwargs)
    # yields output for whatnot.
    yield output

    # keeps track of energy.
    oldenergy = float(structure.energy)
    # gets new structure and iterate.
    structure = crystal.sStructure(output.structure)

  # final calculation.
  vasp.relaxation = "static"
  outdirs.append("%s/final_static" % (outdir))
  output = vasp(structure, outdirs[-1], repat, **kwargs)
  yield output

  # cleanup -- deletes unwanted files from previous output directory
  for dir in outdirs:
    for file in other_repat:
      filename = join(dir, file)
      if exists(filename): remove(filename)


def kpoint_convergence(structure, vasp, outdir="kconv", start=1, steps=None, \
                       offset=(0.5,0.5,0.5), repat=[], tolerance=1e-3, **kwargs):
  """ Performs a convergence test for kpoints using kpoints.Density object.

      This is a generator which yields a output extraction object after each
      vasp calculation:
      >>> for extract in kpoint_convergence(structure, vasp, outdir)
      >>>   print extract.total_energy
      Note that this function will not recompute the result if an appropriate
      output directory is found. 
      @note: This method works only if C{vasp.kpoints = integer(??)} is correct.
      @note: This functor is stateless as long as self and structure can be
             deepcopied correctly.  

      @param structure: A structure to relax.
      @type structure: L{lada.crystal.sStructure} or L{lada.crystal.Structure}.
      @param vasp: The vasp functional.
      @type vasp: L{Vasp}
      @param outdir: Directory where to repatriate files. Default = kconv.
      @type outdir: str
      @param start: Starting density. Default = 1.
      @param steps: End density. Default = 1.
      @param offset: Offset from L{Gamma} of reciprocal mesh. Default = (0.5,0.5,0.5). 
      @param type: 3-tuple.
      @param repat: File to repatriate, other than L{files.minimal}. Default: [].
      @type repat: list or set
      @param tolerance: Total energy convergence criteria. Default: 1e-3. 
      @type tolerance: float
  """
  from copy import deepcopy
  from math import fabs as abs
  from os.path import join, exists
  from os import remove
  from ..kpoints import Density
  from .. import files

  # make this function stateless.
  vasp = deepcopy(vasp)
  repat = set(repat).union(files.minimal)
  tolerance = float(tolerance)
  vasp.relaxation = "static" # what else could it be?
  density = deepcopy(start)

  # keywords arguments cannot include kpoints.
  if "kpoints" in kwargs:
    raise SyntaxError,\
          "Cannot have kpoints as a keyword argument when performing kpoints convergence.\n"

  # we may want to include restart file to speed up calculations. However,
  # these files should be deleted afterwards.
  other_repat = []
  for file in files.restart: 
    if file not in repat: other_repat.append(file)
  repat.union(other_repat)

  # performs initial calculation.
  outdirs = ["%s/density:%i" % (outdir, density)]
  output = vasp(structure, outdirs[-1], repat, kpoints=Density(offset, density), **kwargs)
  # yields output for whatnot
  yield output 
  # makes sure we don't accidentally converge to early.
  oldenergy = -output.total_energy

  tol = tolerance * float( len(structure.atoms) )
  while( abs(oldenergy - output.total_energy) > tol ):
    # restart + new directory
    vasp.indir = outdirs[-1]
    density += steps

    # launch calculations 
    outdirs.append("%s/density:%i" % (outdir, density))
    output = vasp(structure, outdirs[-1], repat, kpoints=Density(offset, density), **kwargs)
    # yields output for whatnot.
    yield output

    # makes sure we don't accidentally converge to early.
    oldenergy = output.total_energy

  # cleanup -- deletes unwanted files from previous output directory
  for dir in outdirs:
    for file in other_repat:
      filename = join(dir, file)
      if exists(filename): remove(filename)
  
def main():
  import os.path
  import shutil
  import subprocess
  from numpy import array as np_array
  from lada import crystal
  from lada.vasp import Vasp, Specie

  structure = crystal.sStructure()
  structure.cell = np_array( [[1.0,0,0],[0,1,0],[0,0,1]] )
  structure.atoms.append( crystal.StrAtom(np_array([0.0,0,0]), "Rb") )
  structure.atoms.append( crystal.StrAtom(np_array([0.5,0.5,0.5]), "K") )
  structure.name = "KRb_nosym"
  structure.scale = 6

  K = Specie( "K", "~/AtomicPotentials/pseudos/K_s" )
  Rb = Specie( "Rb", "~/AtomicPotentials/pseudos/Rb_s" )

  vasp = Vasp(species=(K, Rb))
  vasp.indir = str(structure.name)
  if os.path.exists(vasp.indir):  shutil.rmtree( vasp.indir )

  for extract in relaxation( structure, vasp ):
    structure = extract.structure
    print structure.name, structure.energy, "\n", structure

if __name__ == "__main__":
  main()

    

    
 

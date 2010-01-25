#! /usr/bin/python
""" Method for strain relaxation. """

def relaxation( structure, vasp, outdir, repat = [], tolerance = 1e-3, \
                relaxation="volume ionic cellshape",  noyield = False, **kwargs ):
  """ Performs a vasp relaxation

      This can be called as a generator which yields a output extraction object
      after each vasp calculation:
      >>> for extract in relaxation(structure, vasp)
      >>>   print extract.total_energy
      Or it can be called as a function which yields only a final output extraction object.
      >>> extract = relaxation(structure, vasp, noyield=True).

      Note that this function will not recompute the result if an appropriate
      output directory is found. In other words, this function can be called
      once to perform calculation, and then to perform output extractions.
      >>> # First compute results (if not yet computed).
      >>> relaxation(structure, vasp, "relaxation", noyield = True).
      >>>  ... # do something
      >>> # Then  analyze relaxation steps.
      >>> for extract in relaxation(structure, "relaxation", vasp):
      >>>   print "Where: ", extract.indir
      >>>   print extract.total_energy
      The output extraction object is the output the vasp callable.

      @param structure: A structure to relax.
      @type structure: L{lada.crystal.sStructure} or L{lada.crystal.Structure}.
      @param vasp: The vasp functional.
      @type vasp: L{Vasp}
      @param outdir: Directory where to repatriate files.
      @type outdir: str
      @param repat: File to repatriate, other than L{files.minimal}. Default: [].
      @type repat: list or set
      @param tolerance: Total energy convergence criteria. Default: 1e-3. 
      @type tolerance: float
      @param relaxation: Which degrees of freedom to relax. Default: \"volume
        ionic cellshape\". @see L{incar.Incar.relaxation}.
      @type relaxation: string. 
      @param noyield: whether to use this as a function or a generator. Default: False.
      @type noyield: boolean.
  """
  from copy import deepcopy
  from math import fabs as abs
  from os.path import join, exists, remove
  import files

  # make this function stateless.
  vasp = deepcopy(vasp)
  structure = deepcopy(structure)
  repat = set(repat) + files.minimal
  tolerance = float(tolerance)
  vasp.relaxation = str(relaxation)

  # number of restarts.
  nb_steps = 0 

  # we may want to include restart file to speed up calculations. However,
  # these files should be deleted afterwards.
  other_repat = []
  for file in files.restart: 
    if file not in repat: other_repat.append(file)
  repat += set(other_repat) 

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

    

    
 

#! /usr/bin/env python
from decorator import decorator

@decorator
def count_calls(method, *args, **kwargs):
  """ Adds call counting to a method. """
  result = method(*args, **kwargs)
  args[0]._nbcalls += 1
  return result


def compare_eigenvalues(dft, gw):
  """ Average sum of squares of eigenvalues """
  from numpy import multiply, sum
  mat = dft.eigenvalues - gw.eigenvalues
  return sum( multiply(mat, mat) ) / float(mat.shape[0] * mat.shape[0])

def compare_partial_charges(new, old):
  """ Average sum of squares of partial charges """
  from numpy import multiply, sum
  mat = new.partial_charges - old.partial_charges
  return sum( multiply(mat, mat) ) / float(mat.shape[0] * mat.shape[0])

def compare_pressure(new):
  """ Square of pressure """
  return new.pressure * new.pressure


class Objective(object): 
  """ Objective function to optimize. 

      The vasp object is the one to make actual VASP calls and should be set up
      prior to minimization.
  """
  def __init__(self, vasp, dft, gw, outdir="nlep_fit", comm = None):
    from os import makedirs
    from os.path import exists
    from shutil import rmtree
    from boost.mpi import world
    from lada.crystal import Structure

    self.gw = gw
    self.dft = dft
    self.vasp = vasp
    self.system = Structure(dft.structure)
    self._nbcalls = 0
    self.outdir = outdir
    self.comm = comm
    if self.comm == None: self.comm = world
    if self.comm.rank == 0:
      if exists(self.outdir): rmtree(self.outdir)
      makedirs(self.outdir)
    self.comm.barrier()

  def _get_x0(self):
    """ Returns vector of parameters from L{vasp} attribute. """
    from numpy import array
    result = []
    for specie in self.vasp.species:
      for nlep_params in specie.U:
        if nlep_params["func"] == "nlep": 
          result.append(nlep_params["U"]) 
        elif nlep_params["func"] == "enlep": 
          result.append(nlep_params["U0"]) # first energy
          result.append(nlep_params["U1"]) # second energy
    return array(result, dtype="float64")
  def _set_x0(self, args):
    """ Sets L{vasp} attribute from input vector. """
    from numpy import array
    i = 0
    for specie in self.vasp.species:
      for nlep_params in specie.U:
        if nlep_params["func"] == "nlep": 
          assert args.shape[0] > i+1
          nlep_params["U"] = args[i]
          i += 1
        elif nlep_params["func"] == "enlep": 
          assert args.shape[0] > i+1
          nlep_params["U"] = args[i] # first energy
          nlep_params["J"] = args[i+1] # second energy
          i += 2
    assert self.x == args
  x = property(_get_x0, _set_x0)
  """ Vector of parameters. """

  @count_calls
  def __call__(self, args):
    from os.path import join
    from boost.mpi import world
    from lada.opt.changedir import Changedir
    from lada.vasp import files
    from lada.vasp.extract import Extract
    # transfers parameters to vasp object
    self.x = args
    if world.rank == 0: 
      print join(self.outdir, str(self._nbcalls)), args, self.x
    # performs calculation in new directory
    out = self.vasp\
         (
           self.system,
           outdir = join(self.outdir, str(self._nbcalls)),
           comm = self.comm,
           repat = files.minimal + files.input
         )
    assert out.success, "VASP calculation in %s_%i did not complete." % (self.outdir, self._nbcalls)
    # computes squares
    eigs = compare_eigenvalues(out, self.gw)
    pc = compare_partial_charges(out, self.dft)
    pressure = compare_pressure(out)
    return eigs / 5e0 + pc * 300e0 + pressure / 500e0

  def final(self):
    from os.path import join
    from lada.vasp import Extract
    return Extract(join(self.outdir, str(self._nbcalls))), self.dft, self.gw
      
def main():
  from boost.mpi import world
  from scipy.optimize import fmin as scipy_simplex
  from lada.vasp import Extract, ExtractGW, Vasp, Specie
  from lada.vasp.specie import nlep as nlep_parameters, U as u_parameters
  from lada.vasp.incar import Standard, NBands
  from sys import exit

  indir = "SnO2"
  dft_in = Extract(directory=indir, comm=world)
  dft_in.OUTCAR = "OUTCAR_pbe"
  dft_in.CONTCAR = "POSCAR"
  gw_in = ExtractGW(directory=indir, comm=world)
  gw_in.OUTCAR = "OUTCAR_gw"
  gw_in.CONTCAR = "POSCAR"

  # Creates species with nlep parameters to optimize
  species = Specie\
            (
              "Sn", 
              path="pseudos/Sn", 
              U=[nlep_parameters(type="Dudarev", l=i, U0=0e0) for i in ["s", "p", "d"]]
            ),\
            Specie("O", path="pseudos/O")
  # add U to Sn atoms.
  species[0].U.append( u_parameters(type="Dudarev", U=2e0, l=2) )
  # creates vasp launcher
  vasp = Vasp\
         (
           kpoints    = lambda x: "Automatic generation\n0\ngamma\n6 6 10\n0 0 0",
           precision  = "accurate",
           smearing   = "bloechl",
           ediff      = 1e-5,
           relaxation = "ionic",
           encut      = 1, # uses ENMAX * 1, which is VASP default
           species    = species
         )
  # adds some extra parameters.
  vasp.nbands     = Standard("NBANDS", 64)
  vasp.lorbit     = Standard("LORBIT", 10)
  vasp.npar       = Standard("NPAR", 2)
  vasp.lplane     = Standard("LPLANE", ".TRUE.")
  vasp.addgrid    = Standard("ADDGRID", ".TRUE.")
  del vasp.fftgrid

  # creates objective function.
  objective = Objective(vasp, dft_in, gw_in)

  x0, f0, iter, funcalls, warnflag = scipy_simplex(objective, objective.x, maxfun=5, full_output=1)
  world.barrier()
  if world.rank == 0:
    print "minimum value:", f0
    print "for: ", x0
    print "after %i iterations and %i function calls." % (iter, funcalls)
    print "with warning flag: ", warnflag 

if __name__ == "__main__": main()

    


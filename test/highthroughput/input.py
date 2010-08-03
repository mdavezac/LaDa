vasp = Vasp()
""" VASP functional """
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.set_smearing   = "metal", 0.01
vasp.set_relaxation = "ionic"
vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n30"
# Magmom is a special parameter with specific implementation for point defects.
# We need to go through add_item to reset it as normal key/value parameter.
vasp.add_param      = "magmom", None

# "Al" => specie symbol
# "pseudos/Al" => directory where relevant POTCAR is located
vasp.add_specie = "Al", "pseudos/Al"
vasp.add_specie = "Mg", "pseudos/Mg"
vasp.add_specie =  "O",  "pseudos/O"
vasp.species["Mg"].magnetic = True

first_trial = { "kpoints": "\n0\nAuto\n1", "encut": 0.9 }
""" parameter to override during first relaxation step. """
relaxation_dof = "volume ionic cellshape"
""" Degrees of freedom to relax. """
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)
""" Cell shape relaxation algorithm. """


def volume(structure):
  """ Returns *guessed* scale (eg volume^(0.33)) for a given structure. """
  from numpy.linalg import det
  if "O" in [atom.type for atom in structure.atoms]:    spvol = 8.5**3/4e0
  elif "Se" in [atom.type for atom in structure.atoms]: spvol = 9.5**3/4e0
  elif "Te" in [atom.type for atom in structure.atoms]: spvol = 10.5**3/4e0

  nfu = float(len(structure.atoms)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0) 

def mppalloc(job): 
  """ Returns number of processes for this job. """
  N = len(job.args[0].atoms) # number of atoms.
  if N % 2 == 1: N += 1
  return N  

materials = ["Al2MgO4"]
""" Materials to compute. """
nbantiferro = 3
""" Number of random anti-ferro trials. """


queue     = "regular"
""" PBS queue. """
walltime = "05:45:00"
""" PBS walltime. """

outdir    = "results"
""" Root directory where to store results. """

first_wave_dir = "/scratch/mmdavez/highthroughput"

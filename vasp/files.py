""" Namespace for standard list of files to repatriate. """

POSCAR  = "POSCAR"
""" Name of the input structure file. """
KPOINTS = "KPOINTS"
""" Name of the kpoints file. """
INCAR   = "INCAR"
""" Name of the input parameters file. """
POTCAR  = "POTCAR"
""" Name of the pseudopotential file. """
WAVECAR = "WAVECAR"
""" Name of the wavefunction file. """
CONTCAR = "CONTCAR"
""" Name of the output structure file. """
CHGCAR  = "CHGCAR"
""" Name of the output charge file. """
OSZICAR = "OSZICAR"
""" Name of the energy minimization file. """
STDOUT  = "stdout"
""" Name of the standard output file. """
STDERR  = "stderr"
""" Name of the standard error file. """
EIGENVALUES = "EIGENVALUE"
""" Name of the file with eigenvalues. """
OUTCAR = "OUTCAR"
""" Name of the output file. """

minimal = [OUTCAR, STDOUT, STDERR]
""" The minimum number of files to still have the results of the run.

    minimal = [L{OUTCAR}, L{STDOUT}, L{STDERR}]
"""
output = minimal + [CONTCAR, OSZICAR, EIGENVALUES] 
""" L{Repatriation}.minimal + extra output file. 

    output = minimal + [L{CONTCAR}, L{OSZICAR}, L{EIGENVALUES}] 
"""
input = [INCAR, KPOINTS, POSCAR, POTCAR] 
""" input files. 

    L{restart} = [ L{INCAR}, L{KPOINTS}, L{POSCAR}, L{POTCAR}]
"""
restart = [CHGCAR, WAVECAR]
""" L{Repatriation.safe} + CHGCAR  

    L{restart} = [ L{CHGCAR} + L{WAVECAR} ]
"""
safe = output + input
""" All files to make sense of a run.

    L{safe} = L{Repatriation.output} + L{Repatriation.output} + L{Repatriation.restart} 
"""
all = output + input + restart
""" L{Repatriation.output} + L{Repatriation.output} + L{Repatriation.restart} """

def has_minimal(input):
  """ Checks that the minimum number of files is in repatriation set. """
  for f in minimal:
    if f not in input: return False
  return True


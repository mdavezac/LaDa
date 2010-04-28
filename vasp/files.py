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
EIGENVALUES = "EIGENVAL"
""" Name of the file with eigenvalues. """
OUTCAR = "OUTCAR"
""" Name of the output file. """

minimal = [OUTCAR, CONTCAR, STDOUT, STDERR]
""" The minimum number of files to still have the results of the run.

    minimal = [L{OUTCAR}, L{CONTCAR}, L{STDOUT}, L{STDERR}]
"""
output = minimal + [OSZICAR, EIGENVALUES] 
""" L{minimal} + extra output file. 

    output = L{minimal} + [L{OSZICAR}, L{EIGENVALUES}] 
"""
input = [INCAR, KPOINTS, POSCAR, POTCAR] 
""" input files. 

    L{restart} = [ L{INCAR}, L{KPOINTS}, L{POSCAR}, L{POTCAR}]
"""
restart = [CHGCAR, WAVECAR]
""" L{safe} + L{CHGCAR}

    L{restart} = [ L{CHGCAR} + L{WAVECAR} ]
"""
safe = output + input
""" All files to make sense of a run.

    L{safe} = L{output} + L{input} 
"""
all = output + input + restart
""" L{output} + L{input} + L{restart} """

def has_minimal(files):
  """ Checks that the minimum number of files is in repatriation set. """
  for f in minimal:
    if f not in files: return False
  return True


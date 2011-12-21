""" VASP parameters for lada. """
if "vasp" in globals()["ladamodules"]:
  launch_vasp_as_library = True
  """ Wether to launch vasp as library or program. """
  vasp_program = "vasp"
  """ Path of vasp binary executable (if launching as external program). """

  represent_structure_with_POSCAR = False
  """ If true, then structures are represented using POSCAR format. 

      If False, then uses normal python representation.
  """

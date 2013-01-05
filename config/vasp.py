""" VASP parameters for pylada. """
if "vasp" in globals()["pyladamodules"]:
  vasp_program = "vasp"
  """ Path of vasp binary executable (if launching as external program). """
  vasp_has_nlep = False
  """ Should be set to True if one wants to use NLEP. """


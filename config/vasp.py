""" VASP parameters for pylada. """
if "vasp" in globals()["pyladamodules"]:
  vasp_program = "/projects/nrel/apps/vasp/5.2/bin/vasp"
  """ Path of vasp binary executable (if launching as external program). """
  vasp_has_nlep = False
  """ Should be set to True if one wants to use NLEP. """


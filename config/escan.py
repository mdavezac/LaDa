""" Sets escan parameters. """
if "escan" in globals()["ladamodules"]:
  genpot_library = "libgenpot.so"
  """ Default genpot library. 

      The value for the default can be overriden by ~/.lada in the code below.
  """

  escan_library = "libpescan.so"
  """ Default escan library. 

      The value for the default can be overriden by ~/.lada in the code below.
  """
  launch_escan_as_library = True
  """ Wether to launch escan/genpot as library or program. """
  escan_program = "escanCNL"
  """ Path of escan binary (only needed if launching as external program). """
  genpot_program = "getVLarg"
  """ Path of genpot binary (only needed if launching as external program). """


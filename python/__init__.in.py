""" Lamarck-Darwin extension library."""
import math
import opt
import crystal
import physics
# @do_import_pcm@
@do_import_ce@
@do_import_vasp@
@do_import_escan@
@do_import_vff@
# @do_import_separables@
# @do_import_atompot@
@do_import_enumeration@
@do_import_minimizer@
@do_import_jobs@

version_info = (@LaDa_VERSION_MAJOR@, @LaDa_VERSION_MINOR@)
""" Tuple containing version info. """
version = "%s.%s" % version_info
""" String containing version info. """

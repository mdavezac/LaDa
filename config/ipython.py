""" IPython related configuration. """
auto_import_modules = []
""" Modules to import when starting ipython. """
try_import_matplotlib = False
""" Whether to try and import matplotlib or not. 

    It seems that matplotlib is installed hopper, with the dire consequences
    one expects from Cray.
"""
is_interactive = False
""" Whether we are in an IPython shell or not.

    Tries and converts exceptions to simply printing error messages.
"""
jobdict_glob  = ['*.dict']
""" Globs (unix path wildcards) to identify jobdictionaries. """

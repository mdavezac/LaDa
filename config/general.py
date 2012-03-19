""" Sets general lada parameters. """
jobparams_readonly = False
""" Whether items can be modified in parallel using attribute syntax. """
jobparams_naked_end = True
""" Whether last item is returned as is or wrapped in ForwardingDict. """
jobparams_only_existing = True
""" Whether attributes can be added or only modified. """
unix_re  = True
""" If True, then all regex matching is done using unix-command-line patterns. """
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

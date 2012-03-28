""" Sets general lada parameters. """
jobparams_readonly = False
""" Whether items can be modified in parallel using attribute syntax. """
jobparams_naked_end = True
""" Whether last item is returned as is or wrapped in ForwardingDict. """
jobparams_only_existing = True
""" Whether attributes can be added or only modified. """
unix_re  = True
""" If True, then all regex matching is done using unix-command-line patterns. """
verbose_representation = True
""" Whether functional should be printed verbosely or not. """
ipython_verbose_representation = False
""" When in ipython, should we set :py:data:`verbose_representation` to False. """

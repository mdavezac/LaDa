.. _install_ipython_ug:

Setting up the ipython interface
================================


There are at present three behaviorial variables.

  - :py:data:`pylada.unix_re` is a boolean which controls whether ``jobparams``
    and ``collect`` accept unit-shell-like indices, or whether these are
    regular expressions.
  - :py:data:`pylada.verbose_representation`: Whether objects should be
    transformed to code in a verbose manner, including all attributes with
    their default values, or whether only those attributes which have changed
    from their default values should be explicitely mentionned. When evaluated,
    these two representations should yield the same object. It is recommended
    to keep this value to True, and mess only with the next variable.
  - :py:data:`pylada.ipython_verbose_representation`: If None, this variable is
    ignored. If either True or False, then it
    `:py:data:pylada.verbose_representation` will be set to this value, but only
    within the IPython interface. This means that calculations will represent
    objects in their full glory (eg as pasted at the tail of an OUTCAR), but
    keep the representation short within the interactive interface itself.

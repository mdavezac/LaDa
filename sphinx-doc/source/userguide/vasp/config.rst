.. _vasp_config_ug: 

Configuring LaDa for VASP
=========================

LaDa calls VASP_ as an external program: all it needs is to know how to call it.
The path to the program can be configured in your ~/.lada file by simply adding
the following two lines:

>>> vasp_program = "/path/to/vasp" 
>>> is_vasp_4 = True

:py:data:`~lada.vasp_program` can be absolute path, or simply the name of the
VASP_ binary if it is available in you ``PATH`` environment variable.
:py:data:`~lada.is_vasp_4` should be set to True or False depending on which
version of VASP_ is available. It will prevent some vasp-5 only parameters from
being set and will preferentially write the POSCAR_ in a vasp-5 format.

.. note::

   It is also possible to tell the :py:class:`~functional.Vasp` object to use a
   specific version:
   
   >>> vasp = Vasp(program='/path/to/vasp')
   
   It will apply only to calculations launched with that particular instance.
   If vasp does not have a ``program`` attribute, then it uses the global
   definition.


Since VASP_ is a parallel code, LaDa must also know how to launch MPI binaries.
This is configured for all parallel jobs through :py:data:`~lada.mpirun_exe`:

>>> mpirun_exe = "mpirun -n {n} -npernode {ppn} {program} {cmdline}"

Any string can be entered as long as it possesses the ``{program}`` and
``{cmdline}`` sections. It is a python `format string`_ which is later called
with as user-defined dictionary. However, some options may well be the same for
one run to another. As such a default dictionary can be provided
(:py:data:`~lada.default_comm`).

The last configuration variable is :py:data:`~lada.verbose_representation`. It
controls whether or not the representation/print-out of the functional should
include parameters which have not changed from the default. It is safest to
keep it True.

.. seealso:: :ref:`lada-config`

.. _format string: http://docs.python.org/library/string.html#string-formatting

.. _vasp_running_ug:

Executing a VASP calculation
============================

This is simply the case of calling VASP_ with a structure as input:

>>> vasp(structure, outdir='mycalc')

This will execute VASP_ with the current parameters, in the mycalc directory,
as an mpi process defined by :py:data:`~lada.mpirun_exe` and
:py:data:`lada.default_comm`. The directory can be given using the usual unix
short-cuts and/or shell environment variables.

To specify how the mpi process should go, add a dictionary called ``comm``:

>>> vasp(structure, outdir='~/$WORKDIR/mycalc', comm={'n': 8, 'ppn': 2})

Exactly what this dictionary should contain depends on the specific
supercomputer. The call is formatted by the user-defined
:py:data:`~lada.mpirun_exe`. If the ``comm`` argument is *not* provided, it
defaults to ``None`` and a serial calculation is performed.

Finally, vasp parameters can be modified on a one-off basis:

>>> vasp(structure, outdir='~/mycalc', ispin=1)

.. note::
  
   LaDa will *not* overwrite a successfull calculation, unless specifically
   requested with ``overwrite=True`` in the call, not even one performed
   without the use of LaDa. 

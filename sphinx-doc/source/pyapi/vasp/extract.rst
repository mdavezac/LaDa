Extraction classes
******************
.. module:: lada.vasp.extract

Instances of the extraction classes are returned by calls to the
:py:class:`vasp <lada.vasp.functional.Vasp>` and affiliated methods. They
simply grep the OUTCAR for values of interest, e.g. eigenvalues or vasp
parameters. Indeed, these should contain an ``Extract`` attribute which refers
to a class capable of handling the output of the method or vasp object. They
can be instanciated as follows:

>>> vasp.Extract(outdir)
>>> vasp.epitaxial.Extract(outdir)

Where outdir is the location of the relevant calculation.

The extraction classes are separated into I/O (:py:class:`IOMixin <lada.vasp.extract.mixin.IOMixin>`) and actual
methods to grep the OUTCAR for results (:py:class:`ExtractBase <lada.vasp.extract.base.ExtractBase>`). This setup makes it convenient to change the
kind of object that can be grepped, from the standard file on a hard-disk, to
a file in a database.

:py:class:`Extract` derives from the following classes.

.. toctree::
   :maxdepth: 1

   lada.vasp.extract.base <extractbase>
   lada.vasp.extract.mixin <mixin>


.. autoclass:: Extract
   :members:
   :inherited-members:
  


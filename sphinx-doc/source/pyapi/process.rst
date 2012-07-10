==============
Process Module
==============

.. automodule:: lada.process

  .. moduleauthor:: Mayeul d'Avezac

.. toctree::
   :maxdepth: 1

   Abstract base-class <process/process>
   
   Executing an external program <process/program>

   Executing a generator process <process/iterator>

.. autoclass:: ProcessError
   :show-inheritance:

.. autoclass:: Fail
   :show-inheritance:

.. autoclass:: AlreadyStarted
   :show-inheritance:

.. autoclass:: NotStarted
   :show-inheritance:

.. autofunction:: which

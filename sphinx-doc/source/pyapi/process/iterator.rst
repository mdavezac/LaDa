Executing a generator process
*****************************

.. currentmodule:: lada.process.iterator
.. moduleauthor:: Mayeul d'Avezac
.. autoclass:: IteratorProcess
   :show-inheritance:
   :members:

   .. attribute:: process

      Holds currently running process. 

      This would be the latest process yielded by the input generator
      :py:attr:`functional`.

   .. automethod:: _cleanup

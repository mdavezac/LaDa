lada.vasp.extract.mixin
***********************
.. automodule:: lada.vasp.extract.mixin
.. autoclass:: lada.vasp.extract.mixin.IOMixin
   :members:
   :inherited-members:
   
   .. automethod:: __outcar__()->file object
   .. automethod:: __contcar__()->file object

.. autoclass:: lada.vasp.extract.mixin.OutcarSearchMixin
   :members:
   :inherited-members:

   .. automethod:: _search_OUTCAR(regex, flags=0)->re.match
   .. automethod:: _rsearch_OUTCAR(regex, flags=0)->re.match
   .. automethod:: _find_first_OUTCAR(regex, flags=0)->re.match
   .. automethod:: _find_last_OUTCAR(regex, flags=0)->re.match

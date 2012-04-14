Result mass extraction
**********************

.. automodule:: lada.jobfolder.manipulator
.. moduleauthor:: Mayeul d'Avezac

.. autoclass:: JobParams
   :show-inheritance: 
   :inherited-members:

   .. attribute:: only_existing
      
      If True (default), then only existing parameters can be modified:
      non-existing parameters will not be added.

   .. attribute:: naked_end
      
      If True, then if only one folder contains the requested (sub-)attribute,
      then it is returned as is, rather than wrapped within a
      :py:attr:`~lada.jobfolder.forwarding_dict.ForwardingDict`.

   .. autoattribute:: addattr
   .. autoattribute:: onoff
   .. autoattribute:: extractors
   .. automethod:: __iter__

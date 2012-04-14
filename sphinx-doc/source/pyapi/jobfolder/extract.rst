Extract base classes
********************

.. automodule:: lada.jobfolder.extract
.. moduleauthor:: Mayeul d'Avezac

.. autoclass:: AbstractMassExtract
   :show-inheritance: AbstractMassExtract

   .. automethod:: __init__

   .. autoattribute:: rootpath
   .. autoattribute:: excludes
   .. autoattribute:: view

   .. automethod:: __getitem__
   .. automethod:: __contains__

   .. automethod:: __iter__
   .. automethod:: iteritems
   .. automethod:: itervalues
   .. automethod:: iterkeys
   .. automethod:: items
   .. automethod:: values
   .. automethod:: keys

   .. automethod:: avoid
   .. automethod:: uncache
   .. automethod:: iterfiles

   .. automethod:: shallow_copy

.. autoclass:: AbstractMassExtractDirectories
   :show-inheritance: AbstractMassExtract

   .. automethod:: __init__

   .. automethod:: __is_calc_dir__

Functional class and attributes
*******************************
.. module:: lada.dftcrystal.functional
   :synopsis: Core wrapper for CRYSTAL. 

.. currentmodule:: lada.dftcrystal.functional

.. autoclass:: Functional
   :show-inheritance:
   :members:
   :inherited-members:
   :exclude-members: Extract, iter

   .. autoattribute:: Extract

   .. automethod:: __call__

     This version performs a blocking call to the CRYSTAL_ program. It works as
     follows:

     >>> from lada.dftcrystal import Crystal, Functional
     >>> # create an object
     >>> functional = Functional()
     >>> # initialize it 
     >>> functional.basis['H'] = ...
     ...
     >>> # call the functional on a Crystal structure.
     >>> result = functional(crystal, outdir='ther', comm=comm)
     >>> print result.total_energy
     -666 eV 

   .. automethod:: iter
   
      This version does not perform the actual call to the CRYSTAL_. Rather, it
      is a generator which yields two types of objects:
        
        - :py:attr:`Extract` instances refering to finished CRYSTAL_ calculations
        - :py:mod:`Processes <lada.process>` which allow the user to call the
          CRYSTAL_ code at leisure
     
      The last object yielded should be an extraction object refering to the
      final calculation. In practice, :py:meth:`iter` is used as follows:
   
      .. code-block:: python
       
        for process in functional.iter(crystal, ...):
          # check whether the generator yielded a process or not
          if hasattr(process, 'success'): 
            # This is an extraction object
            # Do something
          else: 
            # This is a process
            # call it somehow
            process.start(comm)
            process.wait()
        print process.some_result # last process should be an extractor

    .. attribute:: dft 

       Interface to the DFT input block of CRYSTAL_.

       It is an instance of :py:class:`hamiltonian.Dft`.
  

.. _CRYSTAL: http://www.crystal.unito.it/

Geometry optimization block
***************************
.. currentmodule:: lada.dftcrystal.optgeom

.. autoclass:: OptGeom
   :show-inheritance:
   :members:
   :inherited-members:
   :exclude-members: fulloptg, cellonly, itatocell, interdun

   .. attribute:: fulloptg

     Optimization of all degrees of freedom: volume, cell-shape, ionic
     positions. 
     
     This is an instance of :py:class:`GeometryOpt`. It excludes
     :py:attr:`cellonly`, :py:attr:`itatocell`, :py:attr:`interdun`.

     >>> functional.optgeom.fulloptg = True
     >>> functional.optgeom.cellonly, functional.optgeom.itatocell, functional.optgeom.interdun
     False, False, False

     It can also be made to optimize at constant volume if :py:attr:`cvolopt` is True.

   .. attribute:: cellonly

     Optimization of cell-shape at constant atomic-position.
     
     This is an instance of :py:class:`GeometryOpt`. It excludes
     :py:attr:`fulloptg`, :py:attr:`itatocell`, :py:attr:`interdun`.

     >>> functional.optgeom.cellonly = True
     >>> functional.optgeom.fulloptg, functional.optgeom.itatocell, functional.optgeom.interdun
     False, False, False

     It can also be made to optimize at constant volume if :py:attr:`cvolopt` is True.

   .. attribute:: itatocell

     Iterative optimization of cell-shape <--> atomic positions. 
     
     This is an instance of :py:class:`GeometryOpt`. It excludes
     :py:attr:`cellonly`, :py:attr:`fulloptg`, :py:attr:`interdun`.

     >>> functional.optgeom.itatocell = True
     >>> functional.optgeom.fulloptg, functional.optgeom.cellonly, functional.optgeom.interdun
     False, False, False

   .. attribute:: interdun

     If True, turns on constrained optimization. See CRYSTAL_ manual.
     
     This is an instance of :py:class:`GeometryOpt`. It excludes
     :py:attr:`cellonly`, :py:attr:`fulloptg`, :py:attr:`interdun`.

     >>> functional.optgeom.itatocell = True
     >>> functional.optgeom.fulloptg, functional.optgeom.cellonly, functional.optgeom.interdun
     False, False, False

   .. attribute:: cvolopt
      
      If True *and* if one of :py:attr:`fulloptg` or :py:attr:`cellonly` is
      True, then performs constant volume optimization.
  
   .. attribute:: bhor
  
      Bhor radius, defined as in CRYSTAL_, 0.5291772083 angstrom.


.. autoclass:: GeometryOpt
   :show-inheritance:
   :members:
   :exclude-members: keyword

   .. attribute:: keyword

      CRYSTAL_ keyword this instance corresponds to.


.. _CRYSTAL: http://www.crystal.unito.it/

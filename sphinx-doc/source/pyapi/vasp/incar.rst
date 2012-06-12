Incar
-----

.. automodule:: lada.vasp.incar
.. moduleauthor:: Mayeul d'Avezac <mayeul.davezac@nrel.gov>

.. autoclass:: Incar
   :show-inheritance:
   :members: 

   Here are all the parameters lada currently has by default.
     :py:attr:`addgrid` :py:attr:`algo` :py:attr:`ediff` :py:attr:`ediffg`
     :py:attr:`encut` :py:attr:`encutgw` :py:attr:`extraelectron`
     :py:attr:`fftgrid` :py:attr:`ispin` :py:attr:`icharg` :py:attr:`istart`
     :py:attr:`isym` :py:attr:`lcharg` :py:attr:`loptics` :py:attr:`lorbit`
     :py:attr:`lmaxmix` :py:attr:`lmaxfockae` :py:attr:`lpead` :py:attr:`lrpa`
     :py:attr:`lsorbit` :py:attr:`lvtot` :py:attr:`lwave` :py:attr:`magmom`
     :py:attr:`nbands` :py:attr:`nomega` :py:attr:`nonscf` :py:attr:`nelm`
     :py:attr:`nelmin` :py:attr:`nelmdl` :py:attr:`npar` :py:attr:`nupdown`
     :py:attr:`precfock` :py:attr:`precision` :py:attr:`relaxation`
     :py:attr:`restart` :py:attr:`symprec` :py:attr:`U_verbosity`
     :py:attr:`system`

   New parameters can be added simply with :py:attr:`add_param`.

   .. include:: incar-attr.rst 

Special Parameters Class Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All these classes are actually defined in ``lada.vasp.incar._params``.

.. autoclass:: Algo
   :show-inheritance:
.. autoclass:: Ediff
   :show-inheritance:
.. autoclass:: Ediffg
   :show-inheritance:
.. autoclass:: Encut
   :show-inheritance:
.. autoclass:: EncutGW
   :show-inheritance:
.. autoclass:: FFTGrid
   :show-inheritance:
.. autoclass:: IniWave
   :show-inheritance:
.. autoclass:: Magmom
   :show-inheritance:
.. autoclass:: NonScf
   :show-inheritance:
.. autoclass:: ExtraElectron
   :show-inheritance:
.. autoclass:: Npar
   :show-inheritance:
.. autoclass:: PartialRestart
   :show-inheritance:
.. autoclass:: PrecFock
   :show-inheritance:
.. autoclass:: Precision
   :show-inheritance:
.. autoclass:: Relaxation
   :show-inheritance:
.. autoclass:: Restart
   :show-inheritance:
.. autoclass:: Smearing
   :show-inheritance:
.. autoclass:: System
   :show-inheritance:
.. autoclass:: UParams
   :show-inheritance:
.. autoclass:: Lsorbit
   :show-inheritance:


.. autoclass:: Choices
    :show-inheritance:
    :members: 

.. autoclass:: Integer
    :show-inheritance:
    :members:

.. autoclass:: Boolean
    :show-inheritance:
    :members:

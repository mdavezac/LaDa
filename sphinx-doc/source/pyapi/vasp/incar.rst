incar
*****

.. automodule:: lada.vasp.incar
.. moduleauthor:: Mayeul d'Avezac <mayeul.davezac@nrel.gov>

.. autoclass:: Incar
   :members: 

   Here are all the parameters lada currently has by default.
     :py:attr:`addgrid` :py:attr:`algo` :py:attr:`ediff` :py:attr:`ediffg`
     :py:attr:`encut` :py:attr:`encutgw` :py:attr:`extraelectron`
     :py:attr:`fftgrid` :py:attr:`ispin` :py:attr:`icharg` :py:attr:`istart`
     :py:attr:`isym` :py:attr:`lcharg` :py:attr:`loptics` :py:attr:`lorbit`
     :py:attr:`lmaxmix` :py:attr:`lmaxfockae` :py:attr:`lpead` :py:attr:`lrpa`
     :py:attr:`lvtot` :py:attr:`lwave` :py:attr:`magmom` :py:attr:`nbands`
     :py:attr:`nomega` :py:attr:`nonscf` :py:attr:`nelm` :py:attr:`nelmin`
     :py:attr:`nelmdl` :py:attr:`npar` :py:attr:`nupdown` :py:attr:`precfock`
     :py:attr:`precision` :py:attr:`relaxation` :py:attr:`restart`
     :py:attr:`symprec` :py:attr:`U_verbosity` :py:attr:`system`

   New parameters can be added simply with :py:attr:`add_param`.

   .. py:attribute:: Incar.algo

      See the documentation for :py:class:`Algo`.

   .. py:attribute:: Incar.addgrid

       Whether to use and additional grid to compute augmentation charges. Defaults to None.

       .. seealso:: `ADDGRID <http://cms.mpi.univie.ac.at/vasp/guide/node142.html>`_

   .. py:attribute:: Incar.ediff

      See the documentation for :py:class:`Ediff`.

   .. py:attribute:: Incar.ediffg

      See the documentation for :py:class:`Ediffg`.

   .. py:attribute:: Incar.encut

      See the documentation for :py:class:`Encut`.

   .. py:attribute:: Incar.encutgw

      See the documentation for :py:class:`EncutGW`.

   .. py:attribute:: Incar.fftgrid

      See the documentation for :py:class:`FFTGrid`.

   .. py:attribute:: Incar.ispin

       Spin-polarized or spin-degenerate calculations.
       Must be 1 or 2. Defaults to 1.

       .. seealso:: `ISPIN <http://cms.mpi.univie.ac.at/vasp/guide/node99.html>`_

   .. py:attribute:: Incar.icharg

       How to intialize the charge density. Defaults to None. 
       It is generally set using :py:attr:`restart <Incar.restart>`.
       Please note that setting :py:attr:`restart <Incar.restart>` to
       something other than None will override what this value.

       .. seealso:: ICHARG_ 
       .. _ICHARG: http://cms.mpi.univie.ac.at/vasp/guide/node102.html

   .. py:attribute:: Incar.istart

       How to intialize the wavefunctions. Defaults to None. 
       It is generally set using :py:attr:`restart <Incar.restart>`.
       Please note that setting :py:attr:`restart <Incar.restart>` to
       something other than None will override what this value.

       .. seealso:: ISTART_
       .. _ISTART: http://cms.mpi.univie.ac.at/vasp/guide/node101.html

   .. py:attribute:: Incar.isym

       Whether and how to use symmetries in the calculation. Defaults to None. 
       You may want to use :py:attr:`symmetries` instead.

       .. seealso:: `ISYM <http://cms.mpi.univie.ac.at/vasp/guide/node115.html>`_

   .. py:attribute:: Incar.lcharg

      Whether to write the charge density to disk. Defaults to True.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LCHARG <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_,
                   :py:Class:`Boolean`

   .. py:attribute:: Incar.lmaxfockae

       Maximum l-quantum number when evaluating on-site terms in the PAW
       method in Hartree fock routines. Defaults to None.

       .. seealso:: `LMAXFOCKAE <http://cms.mpi.univie.ac.at/vasp/vasp/LMAXFOCKAE.html>`_

   .. py:attribute:: Incar.lmaxmix

       Maximum l-quantum number when evaluating on-site terms in the PAW method. Defaults to 4.

       .. seealso:: `LMAXMIX <http://cms.mpi.univie.ac.at/vasp/vasp/PAW_control_tags.html#sec:lmaxmix>`_

   .. py:attribute:: Incar.loptics

      Whether compute frequency dependent dieletric tensor. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

       .. seealso:: `LOPTICS <http://cms.mpi.univie.ac.at/vasp/vasp/LOPTICS_frequency_dependent_dielectric_matrix.html>`_,
                   :py:Class:`Boolean`


   .. py:attribute:: Incar.lorbit

       Whether to write out PROCAR or PROOUT or neither. Defaults to None.

       .. seealso:: `LORBIT <http://cms.mpi.univie.ac.at/vasp/vasp/LORBIT.html>`_
       
   .. py:attribute:: Incar.lpead

      Computes derivative of the wavefunctions with repect to crystal
      momentum using finite difference scheme. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LPEAD <http://cms.mpi.univie.ac.at/vasp/vasp/LPEAD_tag_IPEAD_tag_Derivative_orbitals_w_r_t_k_point.html>`_,
                   :py:Class:`Boolean`

   .. py:attribute:: Incar.lrpa

      Whether to include RPA at the Hartree leval. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LRPA <http://cms.mpi.univie.ac.at/vasp/vasp/LRPA_local_field_effects_on_Hartree_level_RPA.html>`_,
                   :py:Class:`Boolean`

   .. py:attribute:: Incar.lvtot

      Whether to write the charge to disk. Defaults to False.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LVTOT <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_,
                   :py:Class:`Boolean`

   .. py:attribute:: Incar.lwave

      Whether to write wavefunctions to disk. Defaults to False.

      .. seealso:: `LWAVE <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_

   .. py:attribute:: Incar.magmom

      See the documentation for :py:class:`Magmom`.

   .. py:attribute:: Incar.nbands

       Number of bands to consider in the calculation. Defaults to None.

       .. seealso:: `NBANDS <http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html>`_

   .. py:attribute:: Incar.extraelectron

      See the documentation for :py:class:`ExtraElection`.
   
   .. py:attribute:: Incar.nelm

      Maximum number of electronic minimization steps. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELM <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`

   .. py:attribute:: Incar.nelmdl

      Number of non-self-consistent electronic minimization steps at the
      start of the calculation. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELMDL <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`

   .. py:attribute:: Incar.nelmin

      Minimum number of electronic minimization steps. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELMIN <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`


   .. py:attribute:: Incar.nomega

       Number of frequency grid points. Defaults to None.

       .. seealso:: `NOMEGA <http://cms.mpi.univie.ac.at/vasp/vasp/NOMEGA_NOMEGAR_number_frequency_points.html>`_

   .. py:attribute:: Incar.nonscf

      See the documentation for :py:class:`NonScf`.

   .. py:attribute:: Incar.npar

      See the documentation for :py:class:`Npar`.

   .. py:attribute:: Incar.nupdown

       Spin-polarization in number of up eletrons minus number of down electrons. Defaults to None.

       .. seealso:: `NUPDOWN <http://cms.mpi.univie.ac.at/vasp/guide/node122.html>`_

   .. py:attribute:: Incar.precision

      See the documentation for :py:class:`Precision`.

   .. py:attribute:: Incar.precfock

      See the documentation for :py:class:`PrecFock`.

   .. py:attribute:: Incar.relaxation

      See the documentation for :py:class:`Relaxation`.

   .. py:attribute:: Incar.restart

      See the documentation for :py:class:`Restart`.
      .. seealso:: :py:class:`PartialRestart`

   .. py:attribute:: Incar.smearing

      See the documentation for :py:class:`Smearing`.

   .. py:attribute:: Incar.symprec

       Tolerance criteria when determining symmetry operations. Defaults to None. 
       You may want to use :py:attr:`symmetries` instead.

       .. seealso:: `SYMPREC <http://cms.mpi.univie.ac.at/vasp/vasp/ISYM_tag_SYMPREC_tag.html>`_

   .. py:attribute:: Incar.system

      See the documentation for :py:class:`System`.

   .. py:attribute:: Incar.U_verbosity

      See the documentation for :py:class:`UParams`.

Special Parameters Class Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All these classes are actually defined in ``lada.vasp.incar._params``.

.. autoclass:: Algo
.. autoclass:: Ediff
.. autoclass:: Ediffg
.. autoclass:: Encut
.. autoclass:: EncutGW
.. autoclass:: FFTGrid
.. autoclass:: IniWave
.. autoclass:: Magmom
.. autoclass:: NonScf
.. autoclass:: ExtraElectron
.. autoclass:: Npar
.. autoclass:: PartialRestart
.. autoclass:: PrecFock
.. autoclass:: Precision
.. autoclass:: Relaxation
.. autoclass:: Restart
.. autoclass:: Smearing
.. autoclass:: System
.. autoclass:: UParams


.. autoclass:: Choices
    :members: 

.. autoclass:: Integer
    :members:

.. autoclass:: Boolean
    :members:

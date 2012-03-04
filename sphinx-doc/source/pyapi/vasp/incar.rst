incar
*****

.. automodule:: lada.vasp.incar
.. moduleauthor:: Mayeul d'Avezac <mayeul.davezac@nrel.gov>

.. autoclass:: Incar
   :members: 

   Here are all the parameters lada currently has by default.
     :py:attr:`addgrid` :py:attr:`algo` :py:attr:`ediff` :py:attr:`ediffg`
     :py:attr:`encut` :py:attr:`encutgw` :py:attr:`fftgrid` :py:attr:`ibrion`
     :py:attr:`isif` :py:attr:`isigma` :py:attr:`ismear` :py:attr:`ispin`
     :py:attr:`icharg` :py:attr:`istart` :py:attr:`isym` :py:attr:`lcharg`
     :py:attr:`loptics` :py:attr:`lorbit` :py:attr:`lmaxmix`
     :py:attr:`lmaxfockae` :py:attr:`lpead` :py:attr:`lrpa` :py:attr:`lvtot`
     :py:attr:`lwave` :py:attr:`magmom` :py:attr:`nbands` :py:attr:`nomega`
     :py:attr:`nonscf` :py:attr:`nelm` :py:attr:`nelmin` :py:attr:`nelmdl`
     :py:attr:`npar` :py:attr:`nelect` :py:attr:`nsw` :py:attr:`nupdown`
     :py:attr:`potim` :py:attr:`precfock` :py:attr:`precision`
     :py:attr:`restart` :py:attr:`symprec` :py:attr:`U_verbosity`
     :py:attr:`system`

   New parameters can be added simply with :py:attr:`add_param`.

   .. py:attribute:: Incar.ispin

       Spin-polarized or spin-degenerate calculations.
       Must be 1 or 2. Defaults to 1.

       .. seealso:: `ISPIN <http://cms.mpi.univie.ac.at/vasp/guide/node99.html>`_

   .. py:attribute:: Incar.isif

       Which degrees of freedom to relax. I must be 0|1|2|3|4|5|6.
       It can be set more easily *via* :py:attr:`relaxation`. Defaults to 1.

       .. seealso:: `ISIF <http://cms.mpi.univie.ac.at/vasp/guide/node112.html>`_

   .. py:attribute:: Incar.ismear

       Defines the kind of function used for partial occupations. It is
       easier to set using :py:attr:`smearing`. Defaults to None (e.g. to
       VASP default).

       .. seealso:: `ISMEAR <http://cms.mpi.univie.ac.at/vasp/guide/node124.html>`_

   .. py:attribute:: Incar.isigma

       Smearing temperature for partial occupations. It is
       easier to set using :py:attr:`smearing`. Defaults to None (e.g. to
       VASP default).

       .. seealso:: `ISIGMA <http://cms.mpi.univie.ac.at/vasp/guide/node124.html>`_

   .. py:attribute:: Incar.nsw

       Number of ionic steps to perform. Defaults to None.
       Can be set using :py:attr:`relaxation`.

       .. seealso:: `NSW <http://cms.mpi.univie.ac.at/vasp/guide/node108.html>`_

   .. py:attribute:: Incar.ibrion

       Determines the kind of structural relaxation to perform. Defaults to None.
       Can be set using :py:attr:`relaxation`.

       .. seealso:: `IBRION <http://cms.mpi.univie.ac.at/vasp/guide/node110.html>`_

   .. py:attribute:: Incar.potim

       Time step or scaling for forces in some types of molecular dynamics. Defaults to None.
       Can be set using :py:attr:`relaxation`.

       .. seealso:: `POTIM <http://cms.mpi.univie.ac.at/vasp/vasp/POTIM_tag.html>`_

   .. py:attribute:: Incar.nbands

       Number of bands to consider in the calculation. Defaults to None.

       .. seealso:: `NBANDS <http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html>`_

   .. py:attribute:: Incar.lorbit

       Whether to write out PROCAR or PROOUT or neither. Defaults to None.

       .. seealso:: `LORBIT <http://cms.mpi.univie.ac.at/vasp/vasp/LORBIT.html>`_

   .. py:attribute:: Incar.lplan

       Whether real-space data should be distributed plane wise. Defaults to None.

       .. seealso:: `LPLANE <http://cms.mpi.univie.ac.at/vasp/guide/node138.html>`


   .. py:attribute:: Incar.addgrid

       Whether to use and additional grid to compute augmentation charges. Defaults to None.

       .. seealso:: `ADDGRID <http://cms.mpi.univie.ac.at/vasp/guide/node142.html>`_

   .. py:attribute:: Incar.isym

       Whether and how to use symmetries in the calculation. Defaults to None. 
       You may want to use :py:attr:`symmetries` instead.

       .. seealso:: `ISYM <http://cms.mpi.univie.ac.at/vasp/guide/node115.html>`_

   .. py:attribute:: Incar.symprec

       Tolerance criteria when determining symmetry operations. Defaults to None. 
       You may want to use :py:attr:`symmetries` instead.

       .. seealso:: `SYMPREC <http://cms.mpi.univie.ac.at/vasp/vasp/ISYM_tag_SYMPREC_tag.html>`_

   .. py:attribute:: Incar.nupdown

       Spin-polarization in number of up eletrons minus number of down electrons. Defaults to None.

       .. seealso:: `NUPDOWN <http://cms.mpi.univie.ac.at/vasp/guide/node122.html>`_

   .. py:attribute:: Incar.lmaxmix

       Maximum l-quantum number when evaluating on-site terms in the PAW method. Defaults to 4.

       .. seealso:: `LMAXMIX <http://cms.mpi.univie.ac.at/vasp/vasp/PAW_control_tags.html#sec:lmaxmix>`_

   .. py:attribute:: Incar.lmaxfockae

       Maximum l-quantum number when evaluating on-site terms in the PAW
       method in Hartree fock routines. Defaults to None.

       .. seealso:: `LMAXFOCKAE <http://cms.mpi.univie.ac.at/vasp/vasp/LMAXFOCKAE.html>`_

   .. py:attribute:: Incar.nomega

       Number of frequency grid points. Defaults to None.

       .. seealso:: `NOMEGA <http://cms.mpi.univie.ac.at/vasp/vasp/NOMEGA_NOMEGAR_number_frequency_points.html>`_

   .. py:attribute:: Incar.istart

       How to intialize the wavefunctions. Defaults to None. 
       It is generally set using :py:attr:`restart <Incar.restart>`.
       Please note that setting :py:attr:`restart <Incar.restart>` to
       something other than None will override what this value.

       .. seealso:: `ISTART <http://cms.mpi.univie.ac.at/vasp/guide/node101.html>`_

   .. py:attribute:: Incar.icharg

       How to intialize the charge density. Defaults to None. 
       It is generally set using :py:attr:`restart <Incar.restart>`.
       Please note that setting :py:attr:`restart <Incar.restart>` to
       something other than None will override what this value.

       .. seealso:: `ICHARG <http://cms.mpi.univie.ac.at/vasp/guide/node102.html>`_

   .. py:attribute:: Incar.lwave

      Whether to write wavefunctions to disk. Defaults to False.

      .. seealso:: `LWAVE <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_

   .. py:attribute:: Incar.lcharg

      Whether to write the charge density to disk. Defaults to True.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LCHARG <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_,
                   :py:Class:`Boolean`_

   .. py:attribute:: Incar.lvtot

      Whether to write the charge to disk. Defaults to False.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LVTOT <http://cms.mpi.univie.ac.at/vasp/guide/node134.html>`_,
                   :py:Class:`Boolean`_

   .. py:attribute:: Incar.lrpa

      Whether to include RPA at the Hartree leval. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LRPA <http://cms.mpi.univie.ac.at/vasp/vasp/LRPA_local_field_effects_on_Hartree_level_RPA.html>`_,
                   :py:Class:`Boolean`_

   .. py:attribute:: Incar.loptics

      Whether compute frequency dependent dieletric tensor. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

       .. seealso:: `LOPTICS <http://cms.mpi.univie.ac.at/vasp/vasp/LOPTICS_frequency_dependent_dielectric_matrix.html>`_,
                   :py:Class:`Boolean`_

   .. py:attribute:: Incar.lpead

      Computes derivative of the wavefunctions with repect to crystal
      momentum using finite difference scheme. Defaults to None.
      LaDa requires to be True, False, a string representing one the last two,
      or None. 

      .. seealso:: `LPEAD <http://cms.mpi.univie.ac.at/vasp/vasp/LPEAD_tag_IPEAD_tag_Derivative_orbitals_w_r_t_k_point.html>`_,
                   :py:Class:`Boolean`_

   .. py:attribute:: Incar.nelm

      Maximum number of electronic minimization steps. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELM <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`_

   .. py:attribute:: Incar.nelmin

      Minimum number of electronic minimization steps. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELMIN <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`_

   .. py:attribute:: Incar.nelmdl

      Number of non-self-consistent electronic minimization steps at the
      start of the calculation. Defaults to None.
      LaDa requires this argument to be an integer or None.

      .. seealso:: `NELMDL <http://cms.mpi.univie.ac.at/vasp/guide/node104.html>`_,
                   :py:Class:`Integer`_

   .. py:attribute:: Incar.nelect

      See the documentation for :py:class:`NElect`.
   
   .. py:attribute:: Incar.algo

      See the documentation for :py:class:`Algo`.
   
   .. py:attribute:: Incar.precision

      See the documentation for :py:class:`Precision`.

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

   .. py:attribute:: Incar.restart

      See the documentation for :py:class:`Restart`.

   .. py:attribute:: Incar.U_verbosity

      See the documentation for :py:class:`UParams`.

   .. py:attribute:: Incar.magmom

      See the documentation for :py:class:`Magmom`.

   .. py:attribute:: Incar.npar

      See the documentation for :py:class:`Npar`.

   .. py:attribute:: Incar.precfock

      See the documentation for :py:class:`PrecFock`.

   .. py:attribute:: Incar.nonscf

      See the documentation for :py:class:`NonScf`.

   .. py:attribute:: Incar.system

      See the documentation for :py:class:`System`.

Special Parameters Class Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All these classes are actually defined in ``lada.vasp.incar._params``.

.. autoclass:: Restart
.. autoclass:: Algo
.. autoclass:: Encut
.. autoclass:: EncutGW
.. autoclass:: Magmom
.. autoclass:: NElect
.. autoclass:: Precision
.. autoclass:: PrecFock
.. autoclass:: Ediff
.. autoclass:: Ediffg
.. autoclass:: UParams
.. autoclass:: NonScf
.. autoclass:: FFTGrid
.. autoclass:: Npar
.. autoclass:: IniWave
.. autoclass:: System


.. autoclass:: Choices
    :members: 

.. autoclass:: Integer
    :members:

.. autoclass:: Boolean
    :members:

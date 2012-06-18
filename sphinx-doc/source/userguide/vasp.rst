.. _vasp_ug:

Interface to VASP
*****************

LaDa provides an interface wrapping the VASP_ density functional theory code.
This interface manages the input to VASP_, launching the code itself (as an
external program), and retrieving the results as python object. When combined
with the job framework and its attendant ipython_ interface, it becomes a
powerful building block easier to create complex computational schemes.

.. currentmodule:: lada.vasp
    
A fast primer
=============

:py:class:`Vasp <functional.Vasp>` is a python wrapper which hides from the
user all the dreary file manipulation and grepping that working with scientific
codes generally imply. It is initialized as follows:

>>> from lada.vasp import Vasp
>>> vasp = Vasp()
>>> vasp.add_specie = 'Si', '/path/to/directory/Si' 
>>> vasp.add_specie = 'Ge', '/path/to/directory/Ge' 

We will get to the details below, but suppose that we already have some
silicon-germanium structure. Then launching the calculation and retrieving the
gap is as complicated as: 

>>> result = vasp(structure, ispin=2)
>>> print "Success?", result.success
True
>>> print result.eigenvalues
array([-5, -6., -7]) * eV

This will launch vasp in the current directory and wait until the calculation
is finished. We then check whether the calculations completed and print the
eigenvalues. The latter are a numpy array_, signed with the appropriate `units
<qantities>`_.  Since they are a python object, it is easy to create more
complex post-processing.


Configuring LaDa for VASP
=========================

LaDa calls VASP_ as an external program: all it needs is to know how to call it.
The path to the program can be configured in your ~/.lada file by simply adding
the following two lines:

>>> vasp_program = "/path/to/vasp" 
>>> is_vasp_4 = True

:py:data:`vasp_program <lada.vasp_program>` can be absolute path, or simply the
name of the VASP_ binary if it is available in you ``PATH`` environment
variable.  :py:data:`is_vasp_4 <lada.is_vasp_4>` should be set to True or False
depending on which version of VASP_ is available. It will prevent some vasp-5
only parameters from being set and will preferentially write the POSCAR_ in a
vasp-5 format.

.. note:: It is also possible to tell the :py:class:`Vasp <functional.Vasp>`
          object to use a specific version:

          >>> vasp = Vasp(program='/path/to/vasp')

          It will apply only to calculations launched with that particular
          instance. If vasp does not have a ``program`` attribute, then it uses
          the global definition.


Since VASP_ is a parallel code, LaDa must also know how to launch MPI binaries.
This is configured for all parallel jobs through :py:data:`mpirun_exe
<lada.mpirun_exe>`:

>>> mpirun_exe = "mpirun -n {n} -npernode {npernode} {program} {cmdline}"

Any string can be entered as long as it possesses the ``{program}`` and
``{cmdline}`` sections. It is a python `format string`_ which is later called
with as user-defined dictionary. However, some options may well be the same for
one run to another. As such a default dictionary can be provided
(:py:data:`default_comm <lada.default_comm>`).

The last configuration variable is :py:data:`~lada.verbose_representation`. It
controls whether or not the representation/print-out of the functional should
include parameters which have not changed from the default. It is safest to
keep it True.

.. seealso:: :ref:`lada-config`

Creating and Setting up the vasp functional
===========================================

Specifying vasp parameters
--------------------------

As shown in the quick primer above, it is relatively easy to setup a new
:py:class:`~functional.Vasp`. The parameters to the calculation are
generally the VASP_ defaults. They can be changed directly in the constructor:

>>> from lada.vasp import Vasp
>>> vasp = Vasp(ispin=2)

or later, after creation:

>>> vasp.ispin = 1
>>> print vasp.ispin
2

or even right at the moment when executing a calculation:

>>> vasp.encut = 1.2
>>> result = vasp(structure, ispin=2, encut=1)
>>> print vasp.encut
1.2

In the line above, vasp is called with ``ispin=2`` and ``encut=1``. However,
specifying a vasp parameter directly at execution *does* not modify the vasp
object itself, as indicated by the third line. The parameter is only modified
for the duration of the run and not beyond.

At this point, it might have come as a surprise to see ``encut=1``. What units
are we talking about here? A number of parameters have enhanced behaviors with
respect to the original VASP_ code, in order to make it easier to specify
parameters valid for many calculations. :py:attr:`~incar.Incar.encut`
one of those. It can be specified:

  - as floating point smaller than 3, in which case it is a factor of the
    largest `ENMAX` of all the species in the calculation.
  - as floating point larger than 3, in which case it will printed as is in the INCAR.
  - as floating point signed by a quantity:

    >>> from lada.physics import Ry
    >>> vasp.encut = 10*Ry

    In that case, the result will be converted to electron-Volts in the INCAR.

To see the full list of parameters  defined by LaDa, do ``vasp.[TAB]`` on the
ipython_ command line, or go :py:class:`here <lada.vasp.incar.Incar>`. You may
want to checkout :py:attr:`~incar.Incar.relaxation` and
:py:attr:`~incar.Incar.restart`.

The following is a (possibly incomplete) list of VASP parameters with no special behaviors:

.. currentmodule:: lada.vasp.incar

=============================== ================================== ========================== ==========================
:py:attr:`~Incar.addgrid`       :py:attr:`ispin <Incar.ispin>`     :py:attr:`~Incar.istart`   :py:attr:`~Incar.isym`
:py:attr:`~Incar.lmaxfockae`    :py:attr:`lmaxmix <Incar.lmaxmix>` :py:attr:`~Incar.lorbit`   :py:attr:`~Incar.nbands`
:py:attr:`~Incar.nomega`        :py:attr:`nupdown <Incar.nupdown>` :py:attr:`~Incar.symprec` 
=============================== ================================== ========================== ==========================

The following have something special about them:

=================================== ===================================== =================================== =============================== =============================
:py:attr:`~Incar.U_verbosity`       :py:attr:`~Incar.algo`                :py:attr:`~Incar.ediff`             :py:attr:`~Incar.ediffg`        :py:attr:`~Incar.encut`
:py:attr:`~Incar.encutgw`           :py:attr:`~Incar.extraelectron`       :py:attr:`~Incar.fftgrid`           :py:attr:`~Incar.lcharg`        :py:attr:`~Incar.loptics`
:py:attr:`~Incar.lpead`             :py:attr:`~Incar.lrpa`                :py:attr:`~Incar.lsorbit`           :py:attr:`~Incar.lvtot`         :py:attr:`~Incar.lwave`
:py:attr:`~Incar.magmom`            :py:attr:`~Incar.nelm`                :py:attr:`~Incar.nelmdl`            :py:attr:`~Incar.nelmin`        :py:attr:`~Incar.nonscf`
:py:attr:`~Incar.npar`              :py:attr:`~Incar.precfock`            :py:attr:`~Incar.precision`         :py:attr:`~Incar.relaxation`    :py:attr:`~Incar.restart`
:py:attr:`~Incar.smearing`          :py:attr:`~Incar.system`              :py:attr:`~Incar.symmetries`
=================================== ===================================== =================================== =============================== =============================

.. currentmodule:: lada.vasp

Adding parameters LaDa does not know about
------------------------------------------

A large number of parameters control VASP_. LaDa only knows about the most
common. However, it is fairly easy to add parameters which are not there yet.

>>> vasp.add_param = 'encut', 300

This will add a parameter to VASP. It will appear in the incar as "ENCUT =
300", with the tag name in uppercase. The parameter can be later accessed and
modified just as if it where a pre-existing parameter. Indeed, that is how they
are defined in the source code in the first place.

>>> print vasp.encut
300

The attribute ``vasp.encut`` is (in this case) an integer which can be
manipulated as any other python integer. Any python object can be given. It
will appear in the INCAR as it appears above when printed.

If you do *not* want a parameter printed to the INCAR, i.e. you want to use the
VASP_ default for that parameter, the simply set that parameter to None:

>>> vasp.encut = None

.. note:: The examples above specifies a parameter which already 
          exists, :py:attr:`encut <incar.Incar.encut>`.  In this case, when doing
          ``vasp.add_param``,  the parameter is replaced. Its special behavior
          described above is simply gone. If you do not like a special behavior, you
          can always do away with it.

Specifying kpoints
------------------

K-points can be specified in  a variety of ways.

Simply as a KPOINTS_ file in string:

>>> vasp.kpoints = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0" 

As a 4xn matrix where the first three components of each row define a single
kpoint in *cartesian* coordinates, and the last component its multiplicity. In
that case, VASP_ does not reduce the kpoints by symmetry, as explained `here
<KPOINTS>`_.

>>> vasp.kpoints = [ [0, 0, 0, 1], [0.1, 0.1, 0.1, 3] ]


As a callable function which takes the vasp functional and the structure as
input. It must return either a string defining a semantically correct KPOINTS_
file, or an array of kpoints as above (or another similar callable, but why
would you want to do that?).

>>> def some_function(functional, structure):
>>>   .... do something
>>>   return kpoints
>>>
>>> vasp.kpoints = some_function

This function is called before each execution. Hence the kpoint grid can be
tailored to each call.


Specifying the pseudos and Hubbard U parameters
-----------------------------------------------

Pseudo-potentials must be specified explicitely:

>>> vasp.add_specie = 'Tc', '/path/to/pseudo/directory/Tc'

The first item is the name of the atomic specie. It corresponds to the type of
the atom in the structure to compute. The second item is a path to a directory
where the appropriate *unzipped* POTCAR resides. It is possible to specify more
species than exists in a given structure. It will not affect calculations.
For convenience, the path may be given with the usual unix short-cuts and/or with
a shell environment variable.

>>> vasp.add_specie = 'Tc', '~/pseudos/$PAW/Tc'

To specify a Hubbard U parameter, do:

>>> from lada.vasp.specie import U
>>> vasp.add_specie = 'Tc', '~/pseudos/$PAW/Tc', U('dudarev', l=2, U=1.5)

The species can be latter accessed through a dictionary in the vasp object:

>>> vasp.add_specie = 'Tc', '~/pseudos/$PAW/Tc'
>>> print vasp['Tc']
Specie('Tc', '~/pseudos/$PAW/Tc')
>>> vasp['Tc'].moment = 5

At which point, other elemental properties could be added for latter use in a
script.

Executing a VASP calculation
============================

This is simply the case of calling vasp with a structure as input:

>>> vasp(structure, outdir='mycalc')

This will execute VASP_ with the current parameters, in the mycalc directory,
as an mpi process defined by :py:data:`mpirun_exe <lada.mpirun_exe>` and
:py:data:`default_comm <lada.default_comm>`. The directory can be given using
the usual unix short-cuts and/or shell environment variables.

To specify how the mpi process should go, add a dictionary called ``comm``:

>>> vasp(structure, outdir='~/$WORKDIR/mycalc', comm={'n': 8, 'ppn': 2})

Exactly what this dictionary should contain depends on the specific
supercomputer. The call is formatted by the user-defined :py:data:`mpirun_exe
<lada.mpirun_exe>`. If the `comm` argument is *not* provided, it defaults to
:py:data:`default_comm <lada.default_comm>`.

Finally, vasp parameters can be modified on a one-off basis:

>>> vasp(structure, outdir='~/mycalc', ispin=1)

.. note:: LaDa will *not* overwrite a successfull calculation, unless
          specifically requested, not even one performed without the use of
          LaDa. 


Extracting results from a VASP calculation
==========================================

Launching a calculation is not everything. One wants the results. LaDa makes it
easy to interface with those results directly in the python language. Every
functional in LaDa returns an extraction object capable of grepping the
relevant output.

>>> result = vasp(structure)
>>> print result.success
True
>>> print result.eigenvalues * 2
[5, 6, 7] * eV

The above checks that the calculation ran to completion, and then multiplies
the eigenvalues by two. At this point, one could perform any sort of
post-processing, and then automatically launch a subsequent calculation.

.. warning:: Success means the calculation ran to completion, specifically that
             the lines giving the elapsed time exist in the OUTCAR. It does not
             mean that the results are meaningful.
             

The extraction object can be obtained without rerunning VASP_. There are to
ways to go about this. One is, well, to rerun (but without rerunning, if that
makes sense):

>>> # first time, vasp is executed, presumably.
>>> result = vasp(structure, outdir='mycalc', ispin=1)
>>>
>>> # abort if not successful
>>> assert result.success 
>>>
>>> ... do something ...
>>> 
>>> # second time, OUTCAR exists. It is NOT overwritten.
>>> # The extraction object is returned immediately.
>>> result = vasp(structure, outdir='mycalc', ispin=2) 
 
In the example above, vasp is actually launched the first time. However, on the
second pass, an OUTCAR is found. If it is a successful run, then LaDa will
*not* overwrite it. It does not matter whether the structure has changed, or
whether the vasp parameters are different. LaDa will *never* overwrite a
successful run. Not unless specifically requested to. The returned extraction
object corresponds to the OUTCAR. Hence, on the second pass, it is the results
of the first call which are returned. Unless, of course, a successful
calculation already existed there prior to the first run, in which case LaDa
would *never ever* have been so crass as to overwrite it.


The second method is to create an extraction object directly:

>>> from lada.vasp import Extract
>>> result = Extract('/path/to/directory')

In the above, it is expected that the OUTCAR is called OUTCAR. The path can
also be made to a file with a name other than OUTCAR. 

.. note:: An extraction object can be created for any OUTCAR, whether obtained
          *via* LaDa or not. Some information that LaDa automatically appends
          to an OUTCAR may not be obtainable, however.

To find out what LaDa can extract, do ``extract.[TAB]`` in the ipython_
environment. The following is a (possibly incomplete) list of values. Some may
not be available in all calculations, e.g. quasi-particle energies from a DFT
run.

.. currentmodule:: lada.vasp.extract.base

====================================================================================== ====================================================================================== ======================================================================================
:py:attr:`HubbardU_NLEP <ExtractBase.HubbardU_NLEP>`                                   :py:attr:`LDAUType <ExtractBase.LDAUType>`                                             :py:attr:`algo <ExtractBase.algo>`                                                    
:py:attr:`all_total_energies <ExtractBase.all_total_energies>`                         :py:attr:`alphabet <ExtractBase.alphabet>`                                             :py:attr:`cbm <ExtractBase.cbm>`                                                      
:py:attr:`datetime <ExtractBase.datetime>`                                             :py:attr:`density <ExtractBase.density>`                                               :py:attr:`dielectric_constant <ExtractBase.dielectric_constant>`                      
:py:attr:`ediff <ExtractBase.ediff>`                                                   :py:attr:`ediffg <ExtractBase.ediffg>`                                                 :py:attr:`eigenvalues <ExtractBase.eigenvalues>`                                      
:py:attr:`electronic_dielectric_constant <ExtractBase.electronic_dielectric_constant>` :py:attr:`electropot <ExtractBase.electropot>`                                         :py:attr:`encut <ExtractBase.encut>`                                                  
:py:attr:`energies_sigma0 <ExtractBase.energies_sigma0>`                               :py:attr:`energy_sigma0 <ExtractBase.energy_sigma0>`                                   :py:attr:`extraelectron <ExtractBase.extraelectron>`                                  
:py:attr:`fermi0K <ExtractBase.fermi0K>`                                               :py:attr:`fermi_energy <ExtractBase.fermi_energy>`                                     :py:attr:`fft <ExtractBase.fft>`                                                      
:py:attr:`forces <ExtractBase.forces>`                                                 :py:attr:`functional <ExtractBase.functional>`                                         :py:attr:`halfmetallic <ExtractBase.halfmetallic>`                                    
:py:attr:`ialgo <ExtractBase.ialgo>`                                                   :py:attr:`ibrion <ExtractBase.ibrion>`                                                 :py:attr:`icharg <ExtractBase.icharg>`                                                
:py:attr:`initial_structure <ExtractBase.initial_structure>`                           :py:attr:`ionic_charges <ExtractBase.ionic_charges>`                                   :py:attr:`ionic_dielectric_constant <ExtractBase.ionic_dielectric_constant>`          
:py:attr:`is_dft <ExtractBase.is_dft>`                                                 :py:attr:`is_gw <ExtractBase.is_gw>`                                                   :py:attr:`isif <ExtractBase.isif>`                                                    
:py:attr:`ismear <ExtractBase.ismear>`                                                 :py:attr:`ispin <ExtractBase.ispin>`                                                   :py:attr:`istart <ExtractBase.istart>`                                                
:py:attr:`isym <ExtractBase.isym>`                                                     :py:attr:`kpoints <ExtractBase.kpoints>`                                               :py:attr:`lmaxmix <ExtractBase.lmaxmix>`                                              
:py:attr:`lorbit <ExtractBase.lorbit>`                                                 :py:attr:`lvtot <ExtractBase.lvtot>`                                                   :py:attr:`lwave <ExtractBase.lwave>`                                                  
:py:attr:`magnetization <ExtractBase.magnetization>`                                   :py:attr:`moment <ExtractBase.moment>`                                                 :py:attr:`multiplicity <ExtractBase.multiplicity>`                                    
:py:attr:`name <ExtractBase.name>`                                                     :py:attr:`nbands <ExtractBase.nbands>`                                                 :py:attr:`nelmdl <ExtractBase.nelmdl>`                                                
:py:attr:`nelect <ExtractBase.nelect>`                                                 :py:attr:`nelmin <ExtractBase.nelmin>`                                                 :py:attr:`nonscf <ExtractBase.nonscf>`                                                
:py:attr:`nsw <ExtractBase.nsw>`                                                       :py:attr:`nupdown <ExtractBase.nupdown>`                                               :py:attr:`occupations <ExtractBase.occupations>`                                      
:py:attr:`partial_charges <ExtractBase.partial_charges>`                               :py:attr:`potim <ExtractBase.potim>`                                                   :py:attr:`precision <ExtractBase.precision>`                                          
:py:attr:`pressure <ExtractBase.pressure>`                                             :py:attr:`pressures <ExtractBase.pressures>`                                           :py:attr:`pseudopotential <ExtractBase.pseudopotential>`                              
:py:attr:`pulay_pressure <ExtractBase.pulay_pressure>`                                 :py:attr:`qp_eigenvalues <ExtractBase.qp_eigenvalues>`                                 :py:attr:`reciprocal_volume <ExtractBase.reciprocal_volume>`                          
:py:attr:`recommended_fft <ExtractBase.recommended_fft>`                               :py:attr:`relaxation <ExtractBase.relaxation>`                                         :py:attr:`self_energies <ExtractBase.self_energies>`                                  
:py:attr:`sigma <ExtractBase.sigma>`                                                   :py:attr:`smearing <ExtractBase.smearing>`                                             :py:attr:`species <ExtractBase.species>`                                              
:py:attr:`stoichiometry <ExtractBase.stoichiometry>`                                   :py:attr:`stress <ExtractBase.stress>`                                                 :py:attr:`stresses <ExtractBase.stresses>`                                            
:py:attr:`structure <ExtractBase.structure>`                                           :py:attr:`success <ExtractBase.success>`                                               :py:attr:`system <ExtractBase.system>`                                                
:py:attr:`total_energies <ExtractBase.total_energies>`                                 :py:attr:`total_energy <ExtractBase.total_energy>`                                     :py:attr:`valence <ExtractBase.valence>`                                              
:py:attr:`vbm <ExtractBase.vbm>`                                                       :py:attr:`volume <ExtractBase.volume>`                                                 :py:attr:`xc_g0 <ExtractBase.xc_g0>`                                                  

====================================================================================== ====================================================================================== ======================================================================================

.. currentmodule:: lada.vasp

If you know how to use `regular expressions`_, creating a property like those above
is generally fairly simple. Edit the file vasp/extract/base.py, reinstall, and
you're golden. Oh, and send your snippet back this way.

.. _vasp_massextract_ug:

Extracting results from *many* VASP calculation, and plotting stuff
===================================================================

LaDa arranges all calculations within directories, with a single VASP
calculation per sub-directory. It can be expedient to extract simultaneously
all the results contained within a directory and its subdirectories. One
approach is to use :ref:`jobfolders <jobfolder_ug>` and the :ref:`ipython
interface <ipython_ug>`. Failing that, however, it still possible to extract
all the results within a tree of directories simultaneously. When used in
conjunction with plotting software such as matplotlib_, it makes it really easy
to synthesize and understand the results from a set of calculations.
It all comes down to a few simple lines:

>>> from lada.vasp import MassExtract
>>> a = MassExtract('some/path')
>>> a.total_energies
{
   '/this/path/':  array(-666.667) * eV,
   '/that/path/':  array(-999.998) * eV
}

The return is a :py:class:`~lada.jobfolder.forwarding_dict.ForwardingDict`
instance. It is possible to string together attributes to get to those of
interest:

>>> a.structure.scale
{
   '/this/path/':  5.45, 
   '/that/path/':  5.65
}

From there, it is one simple step to plot, say, energies with respect to the
scale (first, run ipython with the ``-pylab`` flag to import matplotlib_
related stuff):

>>> x = array(a.structure.scale.values())
>>> y = array(a.total_energies.values())
>>> plot x, y


:py:class:`~lada.vasp.extract.MassExtract` behaves exactly like the
:ref:`collect <ipython_collect_ug>` object.


Relaxation methods
==================

Relaxing structures generally takes a few actual VASP calculations, since the
FFT and pseudo-potential grids are not adapted as the cell parameters and ionic
positions change. It's brain-dead work best handled automatically. LaDa
currently provides two relaxation methods: :py:class:`relax.Relax` and
:py:class:`relax.Epitaxial`. The former handles general relaxation, including
cellshape, volume, and ionic positions, while the latter performs relaxation on
a virtual substrate (e.g. for coherent Germanium on a (001)@Si substrate, only
the out-of-plane parameter should be relaxed, while the in-plane are fixed to
Si).  The following only describes the first case.

>>> # create the functional
>>> from lada.vasp import Relax
>>> functional = Relax(relaxation='cellshape volume ions', maxcalls=10, keepsteps=True)
>>> functional(structure)

The above creates the functional and launches the calculations. It will first
proceed by relaxing everything, e.g. both cell-shape and ions, if cell-shape
and ionic relaxation are requested. Once convergence is achived, it locks the
cell-shape and relaxes the ions only. Finally, it performs a final static
calculation for maximum accuracy.

The functional is derived from :py:class:`~functional.Vasp`. In practice, this
means that whatever works for :py:class:`~functional.Vasp` works for
:py:class:`~relax.Relax`. However, it does accept a few extra attributes,
described below:

.. glossary::

      first_trial
        A dictionary with parameters which are used only for the very first
        VASP calculation. It can be used to accelerate the first step of the
        relaxation if starting far from the optimum. For instance, it could be
        ``{'encut': 0.8}`` to first converge the structure with a smaller
        cutoff. Defaults to empty dictionary.

      maxcalls
        An interger which denotes the maximum number of calls to VASP before
        aborting. Defaults to 10.

      keepsteps
        If True, intermediate steps are kept. If False, intermediate steps are
        erased.

      relaxation
        Degrees of freedom to relax. It should be either "cellshape" or "ionic"
        or both. Same as for :py:class:`~functional.Vasp`.

      nofail
        If True, will not fail if convergence is not achieved. Just keeps going. 
        Defaults to False.

      convergence
        Convergence criteria. If ``minrelsteps`` is positive, it is only
        checked after ``minrelsteps`` have been performed. Convergence is
        checked according to last VASP run, not from one VASP run to another.
        Eg. If a positive real number, convergence is achieved when the
        difference between the last two total-energies of the current run fall
        below that real number (times structure size), not when the total
        energies of the last two runs fall below that number. Faster, but
        possibly less safe.

        * None: defaults to ``vasp.ediff * 1e1``
        * positive real number: energy convergence criteria in eV per atom. 
        * negative real number: force convergence criteria in eV/angstrom. 
        * callable: Takes an extraction object as input. Should return True if
          convergence is achieved and False otherwise.

      minrelsteps

        Fine tunes how convergence criteria is applied.
        
        * positive: at least ``minrelsteps`` calls to VASP are performed before
          checking for convergence. If ``relaxation`` contains "cellshape",
          then these calls occur during cellshape relaxation. If it does not,
          then the calls occur during the ionic relaxations. The calls do count
          towards ``maxcalls``.
        * negative (default): argument is ignored.


.. note::

   There is also a function :py:func:`relax.relax`, which does the exact same
   thing as the class. Whether to use one or the other is a matter of taste.
   And there's also a generator :py:func:`relax.iter_relax` over actual calls
   to the vasp binary.  Internally, it allows LaDa to launch different
   relaxations in parrallel. In practice, both :py:func:`relax.relax` and
   :py:class:`relax.Relax` make calls to :py:func:`relax.iter_relax`.


Example: Epitaxial relaxation method
====================================

Absent dislocations, a film grown epitaxially on a substrate adopts the
in-plane lattice parameters of that substrate. It is possible to model a
sufficiently thick film as a bulk material on a virtual substrate. The trick is
simply to lock the in-plane parameters while allowing the out-plane direction
to relax.

VASP_ does not allow for this kind of relaxation directly. However, using LaDa
when can simply design a method which will magically perform the relaxation.
Below, a method based on the secant_ method is presented:

  1. Find interval for which stress changes sign.
  2. Bisect interval until required accuracy is obtained.

Admittedly, if one were smart, one could design a faster optimization method.
Or one could use the optimization methods provided by scipy_ with a vasp
object, and be done with it, say this one_.

First, we need a function capable of creating a strain matrix for the
particular kind we have in mind, and apply it to the structure. The strain
matrix can be obtained simply as the outer product of the epitaxial direction with itself.

.. literalinclude:: epirelax.py
   :lines: 1-12

The next step is to create a function which will compute the total energy of a
strained structure while relaxing internal degrees of freedom. 
The vasp call is done in a separate sub-directory.

.. literalinclude:: epirelax.py
   :lines: 14-17

The third component is a function to return the stress component in the
epitaxial direction.

.. literalinclude:: epirelax.py
   :lines: 19-22

Finally, we have all the pieces to put the bisecting algorithm together: 

.. literalinclude:: epirelax.py
   :lines: 24-25, 68-94
   :linenos:

Lines 4 through 16 correspond to (1) above. Starting from the initial input
structure, we change the strain until an interval is found for which the stress
component changes direction. The structure with minimum total energy will
necessarily be found within this interval. Line 14 makes sure the interval is
kept as small as possible.

We then move to the second part (lines 19-25) of the algorithm: bisecting the
interval until the requisite convergence is achieved. Finally, one last static
calculation is performed before returning. That's probably an overkill, but it
should be fairly cheap since we restart from the charge density and
wavefunctions of a previous calculations.

The full listing of the code is given below. There are some extra bits to it,
namely a docstring and some sanity checks on the function's input parameter.

.. literalinclude:: epirelax.py


.. _one: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brent.html
.. _scipy: http://www.scipy.org/
.. _secant: http://en.wikipedia.org/wiki/Secant_method
.. _regular expressions: http://docs.python.org/library/re.html
.. _KPOINTS: http://cms.mpi.univie.ac.at/vasp/guide/node56.html#SECTION00075100000000000000
.. _ENMAX: http://cms.mpi.univie.ac.at/vasp/vasp/ENCUT_tag.html
.. _array: http://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html
.. _quantities: http://packages.python.org/quantities/index.html
.. _ipython: http://ipython.org/
.. _VASP: http://www.vasp.at/
.. _POSCAR: http://cms.mpi.univie.ac.at/vasp/guide/node59.html
.. _format string: http://docs.python.org/library/string.html#string-formatting
.. _matplotlib: http://matplotlib.sourceforge.net/

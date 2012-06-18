.. _install_ug:

Compiling, installing, and setup
********************************

A large subset of the LaDa's functionality needs to be configured for your
specific system, for instance depending on your resource manager, mpi library,
and so forth. In practice, configuration by reading files from one of three places:

 1. Files located in the config sub-directory where LaDa is installed
 2. Files located in one of the directories specified by LADA_CONFIG_DIR
 3. In the user configuration file ~/.lada

Each file is executed and whatever is declared within is placed directly at the
root of the lada package. The files are read in that order. Within a given
directory, files are read alphabetically. Later files will override previous
files, e.g. ~/.lada will override any configuration done previously.

Within an  IPython session, or during an actual calculation, the configuration
variables all reside at the root of the :py:mod:`lada` package.

The following files, located in the ''config'' subdirectory of the source tree,
are examples of different configurations: 

  - cx1.py: PBSpro + intel mpi
  - redmesa.py: SLURM + openmpi


.. toctree::

   Setting up MPI <installation/mpi>

   Setting up the PBS/Slurm ressources <installation/ressources>

   Setting up the ipython interface <installation/ipython>

   Setting up VASP <installation/vasp>

..  Compiling and installing the code
..  =================================
.. 
..  Setting up data directories
..  ===========================
.. 
..  Creating the documentation
..  ==========================

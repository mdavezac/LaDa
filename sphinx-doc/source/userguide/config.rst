Configuring lada
****************

Environment Variables
=====================

.. envvar:: LADA_CONFIG_DIR

   Environment variable specifying the path(s) to the configuration directories.

Configuration variables
=======================

Configuration variables exist in the :py:mod:`lada` module itself. However,
they can be added within separate files. Which files will depend upon the user.

   - Files located in the config sub-directory where lada is installed
   - Files located in one of the directories specified by :envvar:`LADA_CONFIG_DIR`
   - In the user configuration file ~/.lada

Each file is executed and whatever is declared within is placed directly at the
root of the :py:mod:`lada` package. The files are read in that order. Within a
given directory, files are read alphabetically. Later files can override
previous files.

.. currentmodule:: lada

vasp 
----

  These variables are generally declared in config/vasp.py

  .. py:data:: is_vasp_4
     
     If it exists and is True, some vasp parameters will fail if used with
     vasp-5 only options. If it does not exist or is false, then these
     parameters are allowed. 

  .. py:data:: vasp_program

     Path to the vasp executable itself.

===========
Jobs Module
===========

The  :py:mod:`lada.jobfolder` provides tools for high-throughput calculations.
It is centered around an object - the job-folder - which organizes calculations
within a tree of folders, much as one would manually organize calculations
within a tree of directories. Each folder can be executable, e.g.  there is
something to compute there, or non-executable. And each folder can further hold
any number of sub-folders. Furthermore, classes are provided which make it easy
to manipulate the parameters for the calculations in an executable folder, as
well as within all subfolders. Finally, a similar infrastructure is provided to
collect the computational results across all executable sub-folders.


.. toctree::
   :maxdepth: 1

   JobFolder class <jobfolder/jobfolder>

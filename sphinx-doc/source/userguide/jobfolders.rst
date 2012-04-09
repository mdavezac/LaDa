.. _jobfolders_ug: 
.. currentmodule:: lada.jobfolders

Organized high-thoughput calculations: job-folders
**************************************************

LaDa provides tools to organize high-throughput calculations in a systematic
manner. These tools come down to one particular concept: job-folders. A
job-folder is like a directory or folder on a hard-drive. It can contain a
functional and the parameters with which that functional should be called, i.e.
files. And it can contain other job-folders, i.e. sub-directories. Once the
calculations are performed, the architecture of the job-folders is reflected on
disk: the output files from the actual calculations should be found in the
subdirectory corresponding to the sub job-folder.

The following describes how job-folders are created. The fun bits, 
launching jobs, collecting results, manipulating all job-folders
simultaneously, can be found in the next section. Indeed, all of these are
intrinsically linked to the LaDa's IPython interface.


Prep
~~~~
First off, we will need a functional. Rather that use something heavy, like
VASP, we will use a dummy functional which does pretty much nothing... Please
copy the following into a file, any file, which I recommend to call dummy.py.
Putting it into a file is important because we will want python to be able to
refer to it later on.

.. literalinclude:: dummy.py
   :lines: 16-26, 28

This functional takes a few arguments, amongst which an output directory, and
writes a file to disk. That's pretty much it.


Creating and accessing job-folders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Job-folders can be created with two simple lines of codes:

  >>> from lada.jobfolder import JobFolder
  >>> root = JobFolder()

To add further job-folders, one can do:

  >>> jobA = root / 'jobA'
  >>> jobB = root / 'another' / 'jobB'
  >>> jobBprime = root / 'another' / 'jobB' / 'prime'

As you can, see job-folders can be given any structure that on-disk directories
can. What is more, a job-folder can access other job-folders with the same kind
of syntax that one would use (on unices) to access other directories:

  >>> jobA['/'] is root
  True
  >>> jobA['../another/jobB'] is jobB
  True
  >>> jobB['prime'] is jobBprime
  True
  >>> jobBprime['../../'] is not jobB
  True
  >>> root['..']
  KeyError: 'Cannot go below root level.'


Furthermore, job-folders know what they are:

  >>> jobA.name
  '/jobA/'
  
What their parents are:

  >>> jobB.parent.name
  '/another/'

And what the root is:

  >>> jobBprime.root is root
  True
  >>> jobBprime.root.name
  '/'

They also know what they contain:

  >>> 'prime' in jobB
  True
  >>> '/jobA' in jobBprime
  True

Making a job-folder executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The whole point of a job-folder is to create an architecture for calculations.
Each job-folder can contain a single calculation.

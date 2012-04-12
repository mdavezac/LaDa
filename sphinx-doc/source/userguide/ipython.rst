.. currentmodule:: lada
.. _ipython_ug:

IPython high-throughput interface
*********************************

IPython_ is an ingenious combination of a bash-like terminal with a python
shell.  It can be used for both bash related affairs such as copying files
around creating directories, and for actual python programming. In fact, the
two can be combined to create a truly powerfull shell. 

LaDa puts this tool to good use by providing a command-line approach to
manipulate job-folders, launch actual calculations, and collect the result.
When used in conjunction with python plotting libraries, e.g. matplotlib_, it
can provide rapid turnaround from conceptualization to result analysis.

Before moving to the meat, please make sure that LaDa's ipython extension is
running, as :ref:`explained here <install_ip_ug>`.

.. _ipython_prep_ug:

Prep
====

LaDa's IPython interface revolves around :ref:`job-folders <jobfolder_ug>`. 
In order to explore its features, we first need to create job-folders,
preferably some which do not involve heavy calculations. Please copy the
following to a file, such as dummy.py:

.. literalinclude:: dummy.py
   :lines: 1-29, 31-

The above defines three functions: a dummy functional, an extraction function
capable of retrieving the results of the functional from disk, and a function
to create a simple job-folder. In real life, the functional could be a
:py:class:`~vasp.functional.Vasp` object. The extraction function would then be
:py:class:`~vasp.extract.Extract`. And the folder-creation function would depend
on the actual research project.

Manipulating job-folders
========================

Creating a job-folder
~~~~~~~~~~~~~~~~~~~~~

The job-folder could be created as described :ref:`here <jobfolder_ug>`.
However, it is easier -- and safer -- to create a script where a job-folder
creation function resides. If you have performed the step described :ref:`above
<ipython_prep_ug>`, and assuming that the resulting file is called dummy.py,
then a simple job-folder can be created from the ipython interface with:

>>> import dummy
>>> rootfolder = dummy.create_jobs()

.. tip:: 

   We need do ``import`` here because the functional is defined in the script
   itself. However, when using a functional defined somewhere else -- such as
   any LaDa functional -- it is easier for debugging purposes to do:

   >>> run -i dummy.py

   Each time the line above is called, the script is executed anew. When doing
   ``import``, it is executed only the `first time`_. 


Saving and Loading a job-folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At this point we have job-folder stored in memory in a python variable. If you
were to exit ipython, the job-folder would be lost for ever and ever. 

>>> %savefolder dummy.dict rootfolder

The next time ipython is entered, the job-folder can be loaded from disk with: 

>>> %explore dummy.dict

Once a folder has been `explored` from disk, ``savefolder`` can be called
without arguments. 

The percent(%) sign indicates that these commands are ipython
`magic-functions`_. The percent can be obviated using `%automagic`_. To get
more information about what LaDa magic functions do, call them with "--help". 

.. tip::
   
   The current job-folder and the current job-folder path are stored in
   ``lada.interactive.jobfolder`` and ``lada.interactive.jobfolder_path``.
   In practice, accessing those directly is rarely needed.

Listing job-folders
~~~~~~~~~~~~~~~~~~~

Once a job-folder 


.. _IPython: http://ipython.org/
.. _first time: http://ipython.org/ipython-doc/stable/config/extensions/autoreload.html
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _magic-functions: http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions
.. _%automagic: http://ipython.org/ipython-doc/dev/interactive/reference.html#magic-command-system

======================
CRYSTAL wrapper module
======================

.. module:: lada.dftcrystal 
   :synopsis: Wrapper for the CRYSTAL code 

This module creates wrappers around the CRYSTAL_ code. It is main focus are
three kinds of objects: 
  
  - :py:class:`~crystal.Crystal`, :py:class:`~molecule.Molecule`, which define
    the structure to be optimized in a functional manner (as opposed to the
    declarative approach of :py:class:`~lada.crystal.cppwrapper.Structure`)
  - :py:class:`~functional.Functional`, which handles writing the input and
    calling CRYSTAL_ itself
  - :py:class:`~extract.Extract`, which handles grepping values from the output

It strives to reproduce the input style of CRYSTAL_, while still providing a
pythonic interface. As such, the structure CRYSTAL_'s input is mirrored by the
structure of :py:class:`~functional.Functional`'s attributes. 

Content:

.. toctree::

   Functional class and attributes <dftcrystal/functional>

.. currentmodule:: lada.dftcrystal


******
Others
******

.. py:data::  registered

   Map of geometry keywords to their LaDa implementation.


.. py:function:: read(path)->(:py:class:`~crystal.Structure`, :py:class:`~functional.Functional`)

   Reads CRYSTAL_ input from file

   :param path: 
     It can be a file object obtained from open_. Or it could be a string
     defining the path to the input file, or it could a string containing the
     input file itself. The difference in the last two cases resides in whether
     the string has multiple lines or not. Finally, the function is smart
     enough to extract the input from any file, whatever the garbage which
     preceed or follow it. If there are more than one input in the file, only
     the first input is read.

.. _CRYSTAL: http://www.crystal.unito.it/
.. _open: http://docs.python.org/library/functions.html#open



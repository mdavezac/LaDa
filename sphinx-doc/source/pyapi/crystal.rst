The crystal module
==================

.. automodule:: lada.crystal  
.. moduleauthor:: Mayeul d'Avezac <mayeul.davezac@nrel.gov>

.. seealso:: :ref:`crystal_ug`

Classes
-------

Atom
~~~~

.. autoclass:: Atom
  :members: pos, type

  .. automethod:: copy()->Atom
  .. automethod:: to_dict()->dict

Structure
~~~~~~~~~

.. autoclass:: Structure

  .. autoattribute:: volume
  .. autoattribute:: cell
  .. autoattribute:: scale

  .. automethod:: copy()->Structure
  .. automethod:: to_dict()->dict
  .. automethod:: clear()->None
  .. automethod:: insert(index, atom)->None
  .. automethod:: pop(index)->Atom
  .. automethod:: extend(atoms)->None
  .. automethod:: append(atoms)->None
  .. automethod:: transform(matrix)->Structure
  .. automethod:: add_atom(...)->Structure
  .. automethod:: __getitem__()
  .. automethod:: __setitem__()


SmithTransform
~~~~~~~~~~~~~~

.. autoclass:: SmithTransform

  .. autoattribute:: quotient
  .. autoattribute:: transform

  .. automethod:: copy()->SmithTransform
  .. automethod:: indices(position)->numpy.array
  .. automethod:: flatten_indices(indices, index=0)->numpy.array
  .. automethod:: index(position, index=0)->numpy.array


Methods
-------

Vector transformations
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: zero_centered(a, cell, invcell=None)->[numpy.array(..),...]
.. autofunction:: into_voronoi(a, cell, invcell=None)->[numpy.array(..),...]
.. autofunction:: into_cell(a, cell, invcell=None)->[numpy.array(..),...]

Structure transformations
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: supercell(lattice, cell)->Structure
.. autofunction:: primitive(supercell, tolerance=1e-8)->Structure
.. autofunction:: transform(structure, operation)->Structure
.. autofunction:: vasp_ordered

Structure and Geometry Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: is_primitive(supercell, tolerance=1e-8)->Boolean
.. autofunction:: space_group(structure, tolerance=1e-8)->Structure
.. autofunction:: map_sites(mapper, mappee, cmp=None, tolerance=1e-8)->Boolean
.. autofunction:: neighbors(structure, nmax, center, tolerance=1e-8)->[(Atom, numpy.array, float), ...]
.. autofunction:: coordination_shells(structure, nshells, center, tolerance=1e-8, natoms=0)->[ [(Atom, numpy.array, float), ...], ...]
.. autofunction:: periodic_dnc(structure, overlap, mesh=None, nperbox=None, tolerance=1e-8)->[ [(Atom, numpy.array, bool), ...], ...]
.. autofunction:: splitconfigs(structure, center, nmax, configurations=None, tolerance=1e-8)->[ ( (Atom, numpy.array), ...), float), ...]
.. autofunction:: specieset


Iterators
~~~~~~~~~

.. autofunction:: shell_iterator
.. autofunction:: equivalence_iterator
.. autofunction:: layer_iterator


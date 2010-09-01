""" Physical quantities of elements. 

    Atomic quantities can be used as in:

    .. python::
    
      import lada.periodic_table
      periodic_table.Au.atomic_weight
      periodic_table.Gold.atomic_weight

    Available quantities can be found in `Element`. Some of
    these are set to None when either meaningless or not available.
"""
__docformat__ = "restructuredtext en"


from _create_data import *
from _element import Element
import _elements

__all__ = list(_elements.symbols)
__all__.extend(['Element', 'iterate'])
locals().update(_elements.elements)

def iterate():
  """ Iterates through all elements. """
  for name in _elements.symbols:
    yield globals()[name] 


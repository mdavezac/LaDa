//
//  Version: $Id$
//
#ifdef _TYPE_
#  undef _TYPE_
#endif
#ifdef _DIM_
#  undef _DIM_
#endif
#ifdef _CLASSNAME_
#  undef _SHORT_TYPE_ 
#endif
#ifdef _PYTHONNAME_
#  undef _PYTHONNAME_ 
#endif
#ifdef _CONCAT_
#  undef _CONCAT_
#endif
#ifdef _CLASSNAME_
#  undef _CLASSNAME_
#endif

#define _CONCAT_( a, b, c, d ) a ## b ## c ## d

#ifndef _WASINCLUDED_
#  define _WASINCLUDED_ 0
#else
#  undef _WASINCLUDED_
#endif

#if defined(_WASINCLUDED_)

#ifndef _INMODULE_


namespace python
{
  namespace
  {
    std::string print_atom( const ::t_Atom *_at )
    { 
      std::ostringstream sstr;
      _at->print_out( sstr );
      return sstr.str();
    }
  }
}

#else

   class_< t_Atom >( "Atom" )
     .def_readwrite( "pos", &t_Atom::pos )
     .def_readwrite( "site", &t_Atom::site )
     .def_readwrite( "type", &t_Atom::type )
     .def_readwrite( "freeze", &t_Atom::freeze )
     .def( "__str__",  &python::print_atom ) ;

#endif

#include "atom.py.hpp"

#endif 

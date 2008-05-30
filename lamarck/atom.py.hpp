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

namespace Ising_CE
{
}

#else

   class_< Ising_CE::Atom >( "Atom" )
     .def_readwrite( "pos", &Ising_CE::Atom::pos )
     .def_readwrite( "site", &Ising_CE::Atom::site )
     .def_readwrite( "type", &Ising_CE::Atom::type )
     .def_readwrite( "freeze", &Ising_CE::Atom::freeze );

#endif

#include "atom.py.hpp"

#endif 

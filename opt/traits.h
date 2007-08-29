//
//  Version: $Id$
//
#ifndef _TRAITS_H_
#define _TRAITS_H_

#include <list>
#include <vector>

#include <eo/eoPop.h>

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "opt/types.h"
#include "opt/opt_function_base.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Traits 
{
  // As default, T_QUANTITY is expected to be a std style container
  template<class T_QUANTITY>
    struct Quantity
    {
      typedef T_QUANTITY  t_Quantity;  
      typedef typename t_Quantity :: value_type   t_ScalarQuantity;
      typedef std::vector< t_Quantity > t_QuantityGradients;
      bool const static is_vectorial = true;
      bool const static is_scalar = false;
      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned n )
        { return _q[n]; }
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned n )
        { return _q[n]; }
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return _q.size(); }
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::less(_a,_b); }
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::greater(_a,_b); }
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::equal(_a,_b); }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize_container( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_real> 
    {
      typedef types::t_real  t_ScalarQuantity;
      typedef types::t_real  t_Quantity;  
      typedef std::vector< t_Quantity > t_QuantityGradients;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      Quantity() {};
      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return 1; }
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return _a + types::tolerance < _b; }
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return _b + types::tolerance < _a; }
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return std::abs(_b - _a) < types::tolerance; }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_int> 
    {
      typedef types::t_int   t_ScalarQuantity;
      typedef types::t_int   t_Quantity;  
      typedef std::vector< t_Quantity > t_QuantityGradients;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned n ) 
        { return _q; }
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return 1; }
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return _a < _b; }
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return _b < _a; }
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return _b == _a; }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_unsigned> 
    {
      typedef types::t_unsigned   t_ScalarQuantity;
      typedef types::t_unsigned t_Quantity;  
      typedef std::vector< t_Quantity > t_QuantityGradients;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned n ) 
        { return _q; }
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return _a < _b; }
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b )
        { return _b < _a; }
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return _b == _a; }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };

}
#endif

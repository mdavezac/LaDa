//
//  Version: $Id$
//
#ifndef _TRAITS_H_
#define _TRAITS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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
  template< class IS_SCALAR >
   struct Dim { const static bool is_scalar = false;
                const static bool is_vector = true; };
  template<>
   struct Dim<types::t_real> { const static bool is_scalar = true;
                               const static bool is_vector = false; };
  template<>
   struct Dim<types::t_int> { const static bool is_scalar = true;
                              const static bool is_vector = false; };
  template<>
   struct Dim<types::t_unsigned> { const static bool is_scalar = true;
                                   const static bool is_vector = false; };
  template<>
   struct Dim<char> { const static bool is_scalar = true;
                      const static bool is_vector = false; };
  template<>
   struct Dim<bool> { const static bool is_scalar = true;
                      const static bool is_vector = false; };

  template< class T_ARG, bool MAKEVECTOR = true >
   struct MakeVector { typedef std::vector<T_ARG> t_Vector; };
  template< class T_ARG >
   struct MakeVector<T_ARG,false> { typedef T_ARG t_Vector; };

  template< class T_ARG, bool ISVECTOR = Dim<T_ARG>::is_vector >
   struct GetScalar { typedef typename T_ARG::value_type t_Scalar; };
  template< class T_ARG >
   struct GetScalar<T_ARG, false> { typedef T_ARG t_Scalar; };

  template<class T_QUANTITY, bool ISVECTOR = Dim<T_QUANTITY> :: is_vector >
    struct Quantity 
    {
      typedef T_QUANTITY  t_Quantity;  
      typedef typename GetScalar<t_Quantity> :: t_Scalar t_ScalarQuantity;
      typedef Quantity<t_ScalarQuantity>  t_ScalarQuantityTraits;  
      const static bool is_scalar = Dim<t_Quantity> :: is_scalar;
      const static bool is_vector = Dim<t_Quantity> :: is_vector;

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
      static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
      {
        typename t_Quantity :: const_iterator i_scal = _quantity.begin();
        typename t_Quantity :: const_iterator i_end = _quantity.end();
        for(; i_scal != i_end; ++i_scal )
          _stream << *i_scal << " ";
      }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize_container( _q ); }
#endif
    };
  template<class T_QUANTITY >
    struct Quantity <T_QUANTITY, false>
    {
      typedef T_QUANTITY  t_Quantity;  
      typedef typename GetScalar<t_Quantity> :: t_Scalar t_ScalarQuantity;
      typedef Quantity<t_ScalarQuantity>  t_ScalarQuantityTraits;  
      const static bool is_scalar = Dim<t_Quantity> :: is_scalar;
      const static bool is_vector = Dim<t_Quantity> :: is_vector;

      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned n )
        { return _q; }
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return _q.size(); }
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::less(_a,_b); }
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::greater(_a,_b); }
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Quantity<t_ScalarQuantity>::equal(_a,_b); }
      static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
        { _stream << _quantity; }
#ifdef _MPI
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };

  template<class ARGUMENT > class NotRef { typedef ARGUMENT t_notref; };
  template<class ARGUMENT > class NotRef<ARGUMENT& > { typedef ARGUMENT t_notref; };

  template< class T_ARG, class T_RET = T_ARG >
  struct Function
  {
    typedef typename NotRef<T_ARG> :: t_notref  t_Argument;
    typedef typename NotRef<T_RET> :: t_notref  t_Return;
    typedef Quantity<T_ARG> t_ArgTraits;
    typedef Quantity<T_RET> t_RetTraits;
    typedef Quantity< typename MakeVector< t_Return,
                                           Dim<T_ARG>::is_vector > :: t_Vector > t_GradientTraits;
  };

}
#endif

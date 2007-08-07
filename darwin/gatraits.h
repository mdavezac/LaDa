#ifndef _GATRAITS_H_
#define _GATRAITS_H_

#include <list>

#include <eo/eoPop.h>

#include "opt/types.h"
#include "opt/opt_function_base.h"

namespace Traits 
{
  template< class T_INDIVIDUAL,
            class T_OBJECT = typename T_INDIVIDUAL :: t_Object,
            class T_QUANTITY_TRAITS = typename T_INDIVIDUAL :: t_QuantityTraits,
            class T_POPULATION = eoPop< T_INDIVIDUAL >,
            class T_ISLANDS = typename std::list< T_POPULATION > >
  struct Indiv
  {
      typedef T_INDIVIDUAL      t_Individual;
      typedef T_OBJECT          t_Object;
      typedef T_QUANTITY_TRAITS t_QuantityTraits;
      typedef T_POPULATION      t_Population;
      typedef T_ISLANDS         t_Islands;
  };

  template< class T_EVALUATOR,
            class T_INDIVIDUAL = typename T_EVALUATOR :: t_Individual,
            class T_INDIV_TRAITS = typename Traits::Indiv<T_INDIVIDUAL>, 
            class T_QUANTITY_TRAITS = typename T_INDIVIDUAL :: t_QuantityTraits,
            class T_VAFUNCTIONAL = function::Base<types::t_real>, 
            class T_VACONTAINER = typename T_VAFUNCTIONAL :: t_Container,
            class T_VATYPE = typename T_VAFUNCTIONAL :: t_Type >
    struct GA
    {
      typedef T_EVALUATOR t_Evaluator;
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
      typedef T_QUANTITY_TRAITS t_QuantityTraits;
      typedef T_VAFUNCTIONAL  t_VAFunctional;
      typedef T_VACONTAINER t_VAContainer;
      typedef T_VATYPE t_VAType;
    };
//     typedef T_GATRAITS        t_GA_Traits;
//     typedef typename t_GA_Traits :: t_Evaluator        t_Evaluator;
//     typedef typename t_GA_Traits :: t_Individual       t_Individual;
//     typedef typename t_GA_Traits :: t_IndivTraits      t_IndivTraits;
//     typedef typename t_GA_Traits :: t_QuantityTraits   t_QuantityTraits;
//     typedef typename t_GA_Traits :: t_VAFunctional     t_VAFunctional;
//     typedef typename t_GA_Traits :: t_VAContainer      t_VAContainer;
//     typedef typename t_GA_Traits :: t_VAType           t_VAType;
//     typedef typename t_IndivTraits :: t_Population     t_Population;

  template<class T_QUANTITY, class T_CONTAINER = void>
    struct Quantity
    {
      typedef T_CONTAINER  t_Container;  
      typedef T_CONTAINER  t_Quantity;  
      typedef T_CONTAINER& t_Quantity_ref;  
      typedef T_QUANTITY   t_ScalarQuantity;
      typedef T_QUANTITY&  t_ScalarQuantity_ref;
      typedef const T_CONTAINER  const_t_Quantity;  
      typedef const T_CONTAINER& const_t_Quantity_ref;  
      typedef const T_QUANTITY   const_t_ScalarQuantity;
      typedef const T_QUANTITY&  const_t_ScalarQuantity_ref;
      bool const static is_vectorial = true;
      bool const static is_scalar = false;
      Quantity() {};
      t_ScalarQuantity_ref scalar( t_Quantity_ref _q, types::t_unsigned n )
        { return _q[n]; }
      const_t_ScalarQuantity_ref scalar( const_t_Quantity_ref _q, types::t_unsigned n ) const 
        { return _q[n]; }
#ifdef _MPI
      bool broadcast( t_Quantity_ref _q, mpi::BroadCast &_bc )
        { return _bc.serialize_container( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_real> 
    {
      typedef void  t_Container;  
      typedef types::t_real  t_ScalarQuantity;
      typedef types::t_real& t_ScalarQuantity_ref;
      typedef types::t_real  t_Quantity;  
      typedef types::t_real& t_Quantity_ref;  
      typedef const types::t_real  const_t_Quantity;  
      typedef const types::t_real& const_t_Quantity_ref;  
      typedef const types::t_real  const_t_ScalarQuantity;
      typedef const types::t_real& const_t_ScalarQuantity_ref;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      Quantity() {};
      t_ScalarQuantity_ref scalar( t_Quantity_ref _q, types::t_unsigned n )
        { return _q; }
      const_t_ScalarQuantity_ref scalar( const_t_Quantity_ref _q, types::t_unsigned n ) const 
        { return _q; }
#ifdef _MPI
      bool broadcast( t_Quantity_ref _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_int> 
    {
      typedef void  t_Container;  
      typedef types::t_int   t_ScalarQuantity;
      typedef types::t_int&  t_ScalarQuantity_ref;
      typedef types::t_int   t_Quantity;  
      typedef types::t_int&  t_Quantity_ref;  
      typedef const types::t_int  const_t_Quantity;  
      typedef const types::t_int& const_t_Quantity_ref;  
      typedef const types::t_int  const_t_ScalarQuantity;
      typedef const types::t_int& const_t_ScalarQuantity_ref;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      Quantity() {};
      t_ScalarQuantity_ref scalar( t_Quantity_ref _q, types::t_unsigned n )
        { return _q; }
      const_t_ScalarQuantity_ref scalar( const_t_Quantity_ref _q, types::t_unsigned n ) const 
        { return _q; }
#ifdef _MPI
      bool broadcast( t_Quantity_ref _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };
  template<>
    struct Quantity<types::t_unsigned> 
    {
      typedef void  t_Container;  
      typedef types::t_unsigned   t_ScalarQuantity;
      typedef types::t_unsigned&  t_ScalarQuantity_ref;
      typedef types::t_unsigned t_Quantity;  
      typedef types::t_unsigned& t_Quantity_ref;  
      typedef const types::t_unsigned  const_t_Quantity;  
      typedef const types::t_unsigned& const_t_Quantity_ref;  
      typedef const types::t_unsigned  const_t_ScalarQuantity;
      typedef const types::t_unsigned& const_t_ScalarQuantity_ref;
      bool const static is_vectorial = false;
      bool const static is_scalar = true;
      Quantity() {};
      t_ScalarQuantity_ref scalar( t_Quantity_ref _q, types::t_unsigned n )
        { return _q; }
      const_t_ScalarQuantity_ref scalar( const_t_Quantity_ref _q, types::t_unsigned n ) const 
        { return _q; }
#ifdef _MPI
      bool broadcast( t_Quantity_ref _q, mpi::BroadCast &_bc )
        { return _bc.serialize( _q ); }
#endif
    };
}
#endif

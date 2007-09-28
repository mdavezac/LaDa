//
//  Version: $Id$
//
#ifndef _GATRAITS_H_
#define _GATRAITS_H_


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
#include "opt/traits.h"

#include "fitness.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Traits 
{

  template< class T_ARG, class T_RET >
  class VA
  {
      const static bool argvec = Dim<T_ARG> :: is_vector;
      const static bool retvec = Dim<T_RET> :: is_vector;
    public:
      typedef T_ARG                                          t_Container;
      typedef typename GetScalar<t_Container> :: t_Scalar    t_Type;
      typedef function :: Base < t_Type, t_Container >       t_Functional;
                                                                 
    private:                                                     
      typedef typename MakeVector<t_Type, argvec> :: t_Vector  t_A1;

    public:
      typedef typename MakeVector<t_A1, retvec> :: t_Vector  t_QuantityGradients;
  };

  template< class T_OBJECT, class T_CONCENTRATION, class T_FOURIER_RTOK, 
            class T_FOURIER_KTOR = T_FOURIER_RTOK,
            class T_QUANTITY_TRAITS = Traits :: Quantity< typename T_OBJECT :: t_Quantity >,
            class T_VA_TRAITS = Traits::VA<typename T_OBJECT :: t_Container, 
                                           typename T_QUANTITY_TRAITS::t_Quantity >,
            class T_FITNESS = Fitness::Base< T_QUANTITY_TRAITS > >
  struct Indiv
  {
      typedef T_OBJECT                              t_Object;
      typedef T_CONCENTRATION                       t_Concentration;
      typedef T_FOURIER_RTOK                        t_FourierRtoK;
      typedef T_FOURIER_KTOR                        t_FourierKtoR;
      typedef T_QUANTITY_TRAITS                     t_QuantityTraits;
      typedef T_VA_TRAITS                           t_VA_Traits;
      typedef T_FITNESS                             t_Fitness;
      typedef typename t_Fitness :: t_ScalarFitness t_ScalarFitness;
      const static bool is_scalar = t_QuantityTraits :: is_scalar;
      const static bool is_vector = t_QuantityTraits :: is_vector;
  };

  template< class T_EVALUATOR,
            class T_POPULATION = eoPop< typename T_EVALUATOR::t_Individual >,
            class T_ISLANDS = typename std::list< T_POPULATION > >
  struct GA
  {
     typedef T_EVALUATOR   t_Evaluator;
     typedef typename t_Evaluator :: t_Individual   t_Individual;
     typedef T_POPULATION      t_Population;
     typedef T_ISLANDS         t_Islands;
     typedef typename t_Individual :: t_IndivTraits t_IndivTraits;
     typedef typename t_IndivTraits :: t_Object          t_Object;               
     typedef typename t_IndivTraits :: t_Concentration   t_Concentration;
     typedef typename t_IndivTraits :: t_FourierRtoK     t_FourierRtoK;
     typedef typename t_IndivTraits :: t_FourierKtoR     t_FourierKtoR;
     typedef typename t_IndivTraits :: t_QuantityTraits  t_QuantityTraits;
     typedef typename t_IndivTraits :: t_VA_Traits       t_VA_Traits;
     typedef typename t_IndivTraits :: t_Fitness         t_Fitness;
     typedef typename t_Fitness :: t_ScalarFitness       t_ScalarFitness;
     const static bool is_scalar = t_QuantityTraits :: is_scalar;
     const static bool is_vector = t_QuantityTraits :: is_vector;
  };
  

    template< class T_CONTAINER >
      void zero_out( T_CONTAINER &_cont )
      {
        size_t size = _cont.size();
        if ( not size ) return;
        typename T_CONTAINER :: iterator i_first = _cont.begin();
        typename T_CONTAINER :: iterator i_end = _cont.end();
        for(; i_first != i_end; ++i_first)
          zero_out( *i_first );
      }
    template<> inline void zero_out<types::t_real>( types::t_real &_cont )
      { _cont = types::t_real(0);  }
    template<> inline void zero_out<types::t_int>( types::t_int &_cont )
      { _cont = types::t_int(0);  }
    template<> inline void zero_out<types::t_unsigned>( types::t_unsigned &_cont )
      { _cont = types::t_unsigned(0);  }
    template<> inline void zero_out<bool>( bool &_cont )
      { _cont = false;  }
    template<> inline void zero_out<std::string>( std::string &_cont )
      { _cont = "";  }

}
#endif

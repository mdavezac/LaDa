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
      typedef T_ARG                                            t_Container;
      typedef typename opt::GetScalar<t_Container> :: t_Scalar t_Type;
      typedef function :: Base < t_Type, t_Container >         t_Functional;
                                                                 
    private:                                                     
      typedef typename opt::MakeVector<t_Type, argvec> :: t_Vector  t_A1;

    public:
      typedef typename opt::MakeVector<t_A1, retvec> :: t_Vector  t_QuantityGradients;
  };

  template< class T_OBJECT, class T_CONCENTRATION, class T_FOURIER_RTOK, 
            class T_FOURIER_KTOR = T_FOURIER_RTOK,
            class T_QUANTITY_TRAITS = Traits :: Quantity< typename T_OBJECT :: t_Quantity >,
            class T_VA_TRAITS = Traits::VA<typename T_OBJECT :: t_Container, 
                                           typename T_QUANTITY_TRAITS::t_Quantity >,
            class T_FITNESS = typename Fitness::Types< T_QUANTITY_TRAITS > :: Vector >
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
     //! \brief Population of Ref'd individuals
     //! \details Allows to aggregate several populations into one. 
     typedef std::vector< t_Individual* > t_PointerPop; 
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
  
    template< class T_QUANTITY,
              bool is_scalar = Dim< T_QUANTITY > :: is_vector >
      struct SumThem
      {
        explicit
        SumThem( T_QUANTITY &_1, const T_QUANTITY &_2 )
        {
          size_t size1 = _1.size();
          size_t size2 = _2.size();
          if ( size1 == 0 or size2 == 0 ) return;
          if ( size1 != size2 ) return;
        
          typename T_QUANTITY :: iterator i_first = _1.begin();
          typename T_QUANTITY :: const_iterator i_second = _2.begin();
          typename T_QUANTITY :: iterator i_end = _1.end();
          for(; i_first != i_end; ++i_first, ++i_second)
            sum( *i_first, *i_second );
        }
      };
    template< class T_QUANTITY >
      struct SumThem<T_QUANTITY, false>
      {
        explicit
        SumThem( T_QUANTITY &_1, const T_QUANTITY &_2 )
        { _1 = _1 + _2; }
      };
    template<> struct SumThem< std::string >
    {
      explicit
      SumThem( std::string &_1, const std::string &_2 )
        { _1 = _1 + _2; }
    };
    template< class T_QUANTITY >
    void sum( T_QUANTITY& _1, const T_QUANTITY &_2 )
      { SumThem<T_QUANTITY> dummy( _1, _2 ); }

    template< class T_QUANTITY,
              bool is_scalar = Dim< T_QUANTITY > :: is_vector >
      struct ZeroOut
      {
        explicit
        ZeroOut( T_QUANTITY &_cont )
        {
          size_t size = _cont.size();
          if ( not size ) return;
          typename T_QUANTITY :: iterator i_first = _cont.begin();
          typename T_QUANTITY :: iterator i_end = _cont.end();
          for(; i_first != i_end; ++i_first)
            zero_out( *i_first );
        }
      };
    template< class T_QUANTITY >
      struct ZeroOut<T_QUANTITY, false>
      {
        explicit
        ZeroOut( T_QUANTITY &_cont )
          { _cont = T_QUANTITY(0);  }
      };
    template<> struct ZeroOut<std::string>
      {
        explicit
        ZeroOut( std::string &_cont )
          { _cont = "";  }
      };

    template< class T_QUANTITY >
    void zero_out ( T_QUANTITY &_q )
      { ZeroOut<T_QUANTITY> dummy( _q ); }

} // namspace Traits

namespace GA 
{

    template< class T_POPULATION >
      class Aggregate
      {
        public:
          typedef T_POPULATION t_Population; //!< Original population type
        protected:
          //! \cond
          typedef typename t_Population :: value_type t_Individual; 
          typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
          //! \endcond
          //! Aggregate population type
          typedef typename std::vector< t_Individual* > t_Aggregate;

        protected:
          t_Aggregate aggregate;

        public:
          Aggregate() {}; 
          explicit Aggregate( t_Population &_1, t_Population &_2 )
            { ( operator<<( _1 ) ) << _2;  }

          Aggregate& operator<<( t_Population &_pop )
          {
            if( _pop.size() < 1 ) return *this;
            typename t_Population :: iterator i_indiv = _pop.begin();
            typename t_Population :: iterator i_indiv_end = _pop.end();
            for(; i_indiv != i_indiv_end; ++i_indiv )
              aggregate.push_back( &(*i_indiv) );
            return *this;
          }

          operator t_Aggregate& () { return aggregate; }
          operator const t_Aggregate& () const { return aggregate; }
      };
}
#endif

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

#include <mpi/mpi_object.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/traits.h>

#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/traits.h>

#include "fitness.h"


namespace Traits 
{

  //! \ingroup Genetic
  //! \brief Declares input an ouput %types necessary for a Virtual Atom functional.
  //! \param \a T_ARG the type of the argument to the functional
  //! \param \a T_RET the type of the return of the functional
  template< class T_ARG, class T_RET >
  class VA
  {
      //! True if the argument is a vector
      const static bool argvec = Dim<T_ARG> :: is_vector;
      //! true if the return value is a vector
      const static bool retvec = Dim<T_RET> :: is_vector;
    public:
      //! The argument type. See function::Base::t_Container.
      typedef T_ARG                                            t_Container;
      //! The type of the components of the argument. See function::Base::t_Type.
      typedef typename opt::GetScalar<t_Container> :: t_Scalar t_Type;
      //! The type of the functional. 
      typedef function :: Base < t_Type, t_Container >         t_Functional;
                                                                 
    private:                                                     
      //! \cond
      typedef typename opt::MakeVector<t_Type, argvec> :: t_Vector  t_A1;
      //! \endcond

    public:
      //! Type of the gradient of the functional
      typedef typename opt::MakeVector<t_A1, retvec> :: t_Vector  t_QuantityGradients;
  };

  //! \ingroup Genetic
  //! \brief All \e necessary %types for an individual.
  //! \param T_OBJECT declares an object type, e.g. the application specifics
  //!                 which characterize an individual (say a bitstring).
  //! \param T_QUANTITY_TRAITS traits to the quantity used in this
  //!                          optimization. See Traits::Quantity.  
  //! \param T_VA_TRAITS traits for functional minimization with gradients. See Traits::VA.
  //! \param T_FITNESS type of the fitness to use. See Fitness:Types.
  template< class T_OBJECT,
            class T_QUANTITY_TRAITS = Traits :: Quantity< typename T_OBJECT :: t_Quantity >,
            class T_VA_TRAITS = Traits::VA<typename T_OBJECT :: t_Container, 
                                           typename T_QUANTITY_TRAITS::t_Quantity >,
            class T_FITNESS = typename Fitness::Types< T_QUANTITY_TRAITS > :: Vector >
  struct Indiv
  {
      //! Type of the application-specific object 
      typedef T_OBJECT                              t_Object;
      //! \brief Type of the quantity used in search.
      //! \details Determins whether this is a \e scalar or \e vectorial search.
      typedef T_QUANTITY_TRAITS                     t_QuantityTraits;
      //! Defines all %types of a functional, including gradients.
      typedef T_VA_TRAITS                           t_VA_Traits;
      //! Defines the fitness type. 
      typedef T_FITNESS                             t_Fitness;
      //! Defines the scalar fitness. Maybe the same as Indiv::t_Fitness.
      typedef typename t_Fitness :: t_ScalarFitness t_ScalarFitness;
      //! True if the quantity is a scalar.
      const static bool is_scalar = t_QuantityTraits :: is_scalar;
      //! True if the quantity is a vector.
      const static bool is_vector = t_QuantityTraits :: is_vector;
  };

  //! \ingroup Genetic
  //! \brief All \e necessary %types for a \e physical individual.
  //! \param T_OBJECT declares an object type, e.g. the application specifics
  //!                 which characterize an individual (say a bitstring).
  //! \param T_CONCENTRATION is the type to a concentration functor from which
  //!                        the concentration an object and
  //!                        Crystal::Structure can be set and computed.
  //! \param T_FOURIER_RTOK defines the type of a functor which can perform a
  //!                       fourier transform from real space to indirect
  //!                       space. For the arguments to this functor, see
  //!                       Crystal::Fourier.
  //! \param T_FOURIER_RTOK defines the type of a functor which can perform a
  //!                       fourier transform from indirect space to real
  //!                       space. For the arguments to this functor, see
  //!                       Crystal::Fourier.
  //! \param T_QUANTITY_TRAITS traits to the quantity used in this
  //!                          optimization. See Traits::Quantity.  
  //! \param T_VA_TRAITS traits for functional minimization with gradients. See Traits::VA.
  //! \param T_FITNESS type of the fitness to use. See Fitness:Types.
  //! \details To the thypes defined in Traits::Indiv, this class adds a
  //!          concentration functor, a fourier transform from real to indirect
  //!          space, and a fourier transform for indirect to real space. This
  //!          functors are used by the \e physical %GA operators declared in
  //!          gaoperators.h. 
  template< class T_OBJECT, class T_CONCENTRATION, class T_FOURIER_RTOK, 
            class T_FOURIER_KTOR = T_FOURIER_RTOK,
            class T_QUANTITY_TRAITS = Traits :: Quantity< typename T_OBJECT :: t_Quantity >,
            class T_VA_TRAITS = Traits::VA<typename T_OBJECT :: t_Container, 
                                           typename T_QUANTITY_TRAITS::t_Quantity >,
            class T_FITNESS = typename Fitness::Types< T_QUANTITY_TRAITS > :: Vector >
  struct PhysicsIndiv
  {
      //! Type of the application-specific object 
      typedef T_OBJECT                              t_Object;
      //! \brief defines a concentration functor.
      //! \details For the interface expected by this concentration functor,
      //!          see for instance SingleSite::Concentration, or
      //!          Layered::Concentration.
      typedef T_CONCENTRATION                       t_Concentration;
      //! \brief defines a Fourier transform from real to indirect space.
      //! \details It is expected that this fourier transform works well with
      //!          an Crystal::Structure instance. Indeed, for the details of
      //!          the expected interface, see Crystal::Fourier.
      typedef T_FOURIER_RTOK                        t_FourierRtoK;
      //! \brief defines a Fourier transform from indirect to real space.
      //! \details It is expected that this fourier transform works well with
      //!          an Crystal::Structure instance. Indeed, for the details of
      //!          the expected interface, see Crystal::Fourier.
      typedef T_FOURIER_KTOR                        t_FourierKtoR;
      //! \brief Type of the quantity used in search.
      //! \details Determins whether this is a \e scalar or \e vectorial search.
      typedef T_QUANTITY_TRAITS                     t_QuantityTraits;
      //! Defines all %types of a functional, including gradients.
      typedef T_VA_TRAITS                           t_VA_Traits;
      //! Defines the fitness type. 
      typedef T_FITNESS                             t_Fitness;
      //! Defines the scalar fitness. Maybe the same as Indiv::t_Fitness.
      typedef typename t_Fitness :: t_ScalarFitness t_ScalarFitness;
      //! True if the quantity is a scalar.
      const static bool is_scalar = t_QuantityTraits :: is_scalar;
      //! True if the quantity is a vector.
      const static bool is_vector = t_QuantityTraits :: is_vector;
  };

  //! \ingroup Genetic
  //! \brief All relevant %GA %types. 
  //! \details All the %types defined in this structure are necessary for %GA to
  //!          work. Indeed, to compile...
  //! \param T_EVALUATOR type of the interface to application specific behaviors. 
  //! \param T_POPULATION type of the population
  //! \param T_ISLANDS type of the collection of populations.
  template< class T_EVALUATOR,
            class T_POPULATION = eoPop< typename T_EVALUATOR::t_Individual >,
            class T_ISLANDS = typename std::list< T_POPULATION > >
  struct GA
  {
     //! \brief Type of the interface to application specific behaviors. 
     //! \details An instance of this class will enforce \e all application-speciof
     typedef T_EVALUATOR   t_Evaluator;
     //! \brief The type if the individual.
     //! \details Two %types are expected to be declared explicitely in
     //!          t_Evalutor. This is one of them. The other is t_GATraits, eg
     //!          this very structure.
     typedef typename t_Evaluator :: t_Individual   t_Individual;
     //! The type of the container of individuals forming a population
     typedef T_POPULATION      t_Population;
     //! \brief Population of pointers  to individuals
     //! \details Allows to aggregate several populations into one. 
     typedef std::vector< t_Individual* > t_PointerPop; 
     //! A type of container of populations
     typedef T_ISLANDS         t_Islands;
     //! All %types relevant to the %GA individual. See Traits::Indiv.
     typedef typename t_Individual :: t_IndivTraits t_IndivTraits;
     //! Type of the application-specific object. Redeclaration from GA::t_IndivTraits.
     typedef typename t_IndivTraits :: t_Object          t_Object;               
     //! Type of the quantity used in search. Redeclaration from GA::t_IndivTraits.
     typedef typename t_IndivTraits :: t_QuantityTraits  t_QuantityTraits;
     //! \brief Defines all %types of a functional, including gradients. Redeclaration
     //!        from GA::t_IndivTraits.
     typedef typename t_IndivTraits :: t_VA_Traits       t_VA_Traits;
      //! Defines the fitness type. Redeclaration from GA::t_IndivTraits.
     typedef typename t_IndivTraits :: t_Fitness         t_Fitness;
      //! Defines the \e scalar fitness type. Redeclaration from GA::t_IndivTraits.
     typedef typename t_Fitness :: t_ScalarFitness       t_ScalarFitness;
      //! True if the quantity is a scalar. Redeclared from GA::t_IndivTraits.
     const static bool is_scalar = t_QuantityTraits :: is_scalar;
      //! True if the quantity is a vector. Redeclared from GA::t_IndivTraits.
     const static bool is_vector = t_QuantityTraits :: is_vector;
  };
  
  //! \ingroup Genetic
  //! \brief Sums two objects, wether scalars or vectors, or even strings
  //! \details More specifically, this functor takes two arguments \a A and \a
  //!           B, and sums B into A. The functor is declared as a constructor,
  //!           so you can call it like you would any other %function. The
  //!           default implementation is for vectors.
  //!           The %function sum() is a helper %function for this class.
  template< class T_QUANTITY,
            bool is_scalar = Dim< T_QUANTITY > :: is_vector >
    struct SumThem
    {
      //! The functor/Constructor thingee
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
  //! \ingroup Genetic
  //! \brief Sums two objects. Scalar flavor.
  template< class T_QUANTITY >
    struct SumThem<T_QUANTITY, false>
    {
      //! The functor/Constructor thingee
      explicit
      SumThem( T_QUANTITY &_1, const T_QUANTITY &_2 )
      { _1 = _1 + _2; }
    };
  //! \ingroup Genetic
  //! \brief Sums two objects. std::string flavor.
  template<> struct SumThem< std::string >
  {
    //! The functor/Constructor thingee
    explicit
    SumThem( std::string &_1, const std::string &_2 )
      { _1 = _1 + _2; }
  };
  //! \ingroup Genetic
  //! Helper %function to SumThem.
  template< class T_QUANTITY >
  void sum( T_QUANTITY& _1, const T_QUANTITY &_2 )
    { SumThem<T_QUANTITY> dummy( _1, _2 ); }

  //! \ingroup Genetic
  //! \brief Initializes scalars and vectors and strings to zero.
  //! \details More specically, numeric scalars are set to zero, std::string
  //!         instances are set to empty, and the components of a vector are
  //!         set to zero (as defined by the type of the component). The
  //!         functor is declared as a constructor, so you can call it like you
  //!         would any other %function. The default implementation is for
  //!         vectors. The %function zero_out() is a helper %function for this class.
  template< class T_QUANTITY,
            bool is_scalar = Dim< T_QUANTITY > :: is_vector >
    struct ZeroOut
    {
      //! The functor/Constructor thingee
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
  //! \ingroup Genetic
  //! \brief Initializes scalars to zero.
  template< class T_QUANTITY >
    struct ZeroOut<T_QUANTITY, false>
    {
      //! The functor/Constructor thingee
      explicit
      ZeroOut( T_QUANTITY &_cont )
        { _cont = T_QUANTITY(0);  }
    };
  //! \ingroup Genetic
  //! \brief Initializes std::string instance to an empty string.
  template<> struct ZeroOut<std::string>
    {
      //! The functor/Constructor thingee
      explicit
      ZeroOut( std::string &_cont )
        { _cont = "";  }
    };

  //! \ingroup Genetic
  //! \brief Helper %function for ZeroOut functor.
  template< class T_QUANTITY >
  void zero_out ( T_QUANTITY &_q )
    { ZeroOut<T_QUANTITY> dummy( _q ); }

} // namspace Traits

namespace GA 
{
  
    //! \brief creates an aggregate population of pointers to individuals from
    //!        two populations of individuals.
    //! \details This functor can be used quite simply in the following way.
    //! \code
    //    Aggregate<t_Population> aggregate.
    //    aggregate << population1 << population2 << .... << populationN;
    //    std::vector< t_Individual* >& newpop = aggregate.
    //    // do whatever you want with newpop.
    //! \endcode
    //! Note that Aggregate can always be used in place of a reference to
    //! Aggregate::t_Aggregate, so you need not define \a newpop explicitely as
    //! above.
    //! 
    //! You can use Aggregate in conjunction with Modifier::innermost and
    //! Modifier::const_innermost, so that one and the same code can deal with a
    //! collection of pointers, as well as a collection of instances or
    //! references. See for instance the implementation of Scaling::Pareto.
    template< class T_POPULATION >
      class Aggregate
      {
        public:
          typedef T_POPULATION t_Population; //!< Original population type
        protected:
          //! The type of the individual
          typedef typename t_Population :: value_type t_Individual; 
          //! All %%types relevant to the individuals.
          typedef typename t_Individual :: t_IndivTraits t_IndivTraits; 
          //! Aggregate population type
          typedef typename std::vector< t_Individual* > t_Aggregate;

        protected:
           //! \brief The aggregated population.
           //! \details The aggregate population is made of pointers to the
           //!          original instances of the individuals in the original
           //!          population.
          t_Aggregate aggregate; 

        public:
          //! Constructor.
          Aggregate() {}; 
          //! Explicit constructor aggregating two populations.
          explicit Aggregate( t_Population &_1, t_Population &_2 )
            { ( operator<<( _1 ) ) << _2;  }

          //! \brief Appends pointers to the components of \a  _pop at the end
          //!        of Aggregate :: aggregate. 
          Aggregate& operator<<( t_Population &_pop )
          {
            if( _pop.size() < 1 ) return *this;
            typename t_Population :: iterator i_indiv = _pop.begin();
            typename t_Population :: iterator i_indiv_end = _pop.end();
            for(; i_indiv != i_indiv_end; ++i_indiv )
              aggregate.push_back( &(*i_indiv) );
            return *this;
          }

          //! Returns a reference to Aggregate::aggregate.
          operator t_Aggregate& () { return aggregate; }
          //! Returns a constant reference to Aggregate::aggregate.
          operator const t_Aggregate& () const { return aggregate; }
      };
}
#endif

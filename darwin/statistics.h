//
//  Version: $Id$
//
#ifndef _DARWIN_STATISTICS_H_
#define _DARWIN_STATISTICS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoPop.h>
#include <eo/eoGenOp.h>

#include <opt/types.h>
#include <print/xmg.h>

namespace GA
{
  /** \brief Counts the number of non-identical individuals in the population
  * \details What is meant by "non-identical" will depend on the implementation of
  *          t_Individual::operator==() */
  template< class T_GATRAITS>
  class TrueCensus : public eoStatBase< typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< Contains all %GA types
    protected:
      typedef typename t_GATraits::t_Population t_Population; //!< Population type
      typedef typename t_GATraits::t_Individual t_Individual; //!< Individual type
      typedef typename t_GATraits::t_Object t_Object; //!< Object type
      
    public:
      TrueCensus () {} //!< Constructor
      virtual ~TrueCensus() {} //!< Destructor
      //! Name of the class, EO requirement
      virtual std::string className(void) const { return "GA::TrueCensus"; }
      //! \brief Functor which counts the number of unique individuals in a population \a _pop
      //! \details Results are printed out on Print::xmg as a comment
      virtual void operator()( const t_Population &_pop );
  };

  //! \brief Prints average \e scalar fitness over the current population 
  //! \details Repeat. This is the average of the \e scalar fitness, not of the \e
  //!          vectorial fitness. 
  template< class T_GATRAITS>
  class AverageFitness : public eoStatBase< typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< Contains all %GA types
    protected:
      typedef typename t_GATraits::t_Population t_Population; //!< Population type
      typedef typename t_GATraits::t_Individual t_Individual; //!< Individual type
      typedef typename t_GATraits::t_Object t_Object; //!< Object type
      
    public:
      AverageFitness () {} //!< Constructor
      virtual ~AverageFitness() {} //!< Destructor
      //! Name of the class, EO requirement
      virtual std::string className(void) const { return "GA::AverageFitness"; }
      //! Computes the average \e scalar fitness over the currrent population.
      virtual void operator()( const t_Population &_pop );
  };
  

  //! \brief Prints average quantity over the current population 
  //! \details Repeat. This is the average of the \e scalar quantities in the
  //!          case of single-objective %GA, and of vectorial quantities in the
  //!          case of multi-objective %GA.
  template< class T_GATRAITS,
            bool is_scalar = T_GATRAITS::t_QuantityTraits::is_scalar >
  class AverageQuantities : public eoStatBase< typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< Contains all %GA types
    protected:
      typedef typename t_GATraits::t_Population t_Population; //!< Population type
      typedef typename t_GATraits::t_Individual t_Individual; //!< Individual type
      typedef typename t_GATraits::t_Object t_Object; //!< Object type
      
    public:
      AverageQuantities () {} //!< Constructor
      virtual ~AverageQuantities() {} //!< Destructor
      //! Name of the class, EO requirement
      virtual std::string className(void) const { return "GA::AverageQuantities"; }
      //! Computes the average quantity (-ies) over the currrent population.
      virtual void operator()( const t_Population &_pop );
  };
  //! \brief Prints average quantity over the current population. Scalar
  //!        quantities flavor.
  template< class T_GATRAITS >
  class AverageQuantities<T_GATRAITS, true> :
    public eoStatBase< typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< Contains all %GA types
    protected:
      typedef typename t_GATraits::t_Population t_Population; //!< Population type
      typedef typename t_GATraits::t_Individual t_Individual; //!< Individual type
      typedef typename t_GATraits::t_Object t_Object; //!< Object type
      
    public:
      AverageQuantities () {} //!< Constructor
      virtual ~AverageQuantities() {} //!< Destructor
      //! Name of the class, EO requirement
      virtual std::string className(void) const { return "GA::AverageQuantities"; }
      //! Computes the average quantity (-ies) over the currrent population.
      virtual void operator()( const t_Population &_pop );
  };
  


  template< class T_GATRAITS>
  inline void TrueCensus<T_GATRAITS> :: operator()( const t_Population &_pop )
  {
    typename t_Population :: const_iterator i_begin = _pop.begin();
    typename t_Population :: const_iterator i_indiv1 = i_begin;
    typename t_Population :: const_iterator i_end = _pop.end();
    typename t_Population :: const_iterator i_indiv2;
    types::t_unsigned N = _pop.size();
    
    for( ; i_indiv1 != i_end; ++i_indiv1 )
    {
      i_indiv2 = i_indiv1; ++i_indiv2;
      if ( i_end != std::find( i_indiv2, i_end, *i_indiv1 ) )
        --N;
    }

    // prints stuff out
    Print::xmg << Print::Xmg::comment <<  "True Census: "
               << N << " / " << _pop.size() << Print::endl; 
  }



  template< class T_GATRAITS>
  inline void AverageFitness<T_GATRAITS> :: operator()( const t_Population &_pop )
  {
    typename t_Population :: const_iterator i_indiv = _pop.begin();
    typename t_Population :: const_iterator i_end = _pop.end();
    types::t_unsigned N = _pop.size();
    
    typedef typename t_GATraits :: t_QuantityTraits :: t_ScalarQuantity t_S;
    t_S average(0);
    for( ; i_indiv != i_end; ++i_indiv )
      average += (const t_S&) i_indiv->fitness();
    
    average /= (t_S) _pop.size();
    Print::out << "Average Fitness: " << average << "\n\n";
  }


  template< class T_GATRAITS, bool is_scalar>
  inline void AverageQuantities<T_GATRAITS, is_scalar> ::
     operator()( const t_Population &_pop )
  {
    typename t_Population :: const_iterator i_indiv = _pop.begin();
    typename t_Population :: const_iterator i_indiv_end = _pop.end();
    
    typedef typename t_GATraits :: t_QuantityTraits :: t_Quantity t_V;
    typedef typename t_GATraits :: t_QuantityTraits :: t_ScalarQuantity t_S;
    t_V average;
    average.resize( i_indiv->const_quantities().size() );
    Traits::zero_out( average );
    for( ; i_indiv != i_indiv_end; ++i_indiv )
      Traits::sum( average,  i_indiv->const_quantities() );
    
    std::transform( average.begin(), average.end(), average.begin(),
                    std::bind2nd( std::divides<t_S>(), (t_S) _pop.size() ) );

    Print::out << "Average Quantity: "
               <<  Traits::Quantity<t_V>::print( average )
               << "\n\n";
  }
  template< class T_GATRAITS >
  inline void AverageQuantities<T_GATRAITS, true> ::
     operator()( const t_Population &_pop )
  {
    typename t_Population :: const_iterator i_indiv = _pop.begin();
    typename t_Population :: const_iterator i_indiv_end = _pop.end();
    
    typedef typename t_GATraits :: t_QuantityTraits :: t_Quantity t_V;
    types::t_real average(0);
    for( ; i_indiv != i_indiv_end; ++i_indiv )
      average += (t_V) i_indiv->const_quantities();
    
    average /= ( types::t_real ) _pop.size();
    Print::out << "Average Quantity: "
               << average
               << "\n\n";
  }
}
  
#endif

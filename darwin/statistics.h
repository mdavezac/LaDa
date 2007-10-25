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

#include "opt/types.h"
#include "print/xmg.h"

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
      virtual std::string className(void) const { return "GA::TrueCensus"; }
      //! \brief Functor which counts the number of unique individuals in a population \a _pop
      //! \details Results are printed out on Print::xmg as a comment
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
      if ( i_end != std::find_if( i_indiv2, i_end,
                                  std::bind1st(std::equal_to<t_Object>(), *i_indiv1) ) )
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
}
  
#endif

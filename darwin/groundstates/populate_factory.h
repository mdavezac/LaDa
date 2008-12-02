#ifndef _LADA_GA_GROUNDSTATES_POPULATEFACTORY_H_
#define _LADA_GA_GROUNDSTATES_POPULATEFACTORY_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <string>

#include "../operators/populate.h"

namespace LaDa
{
  namespace GA
  {
    namespace GroundStates
    {
      //! Helps bind a populate functor to an attribute factory.
      template< class T_EVALUATOR >
        void populate_factory( const std::string& _string,
              boost::function<void(typename T_EVALUATOR::t_GATraits::t_Population&, 
                                   size_t)>& _func,
              T_EVALUATOR& _evaluator,
              const Taboo::Container<typename T_EVALUATOR::t_GATraits::t_Individual>& _taboo )
        {
          namespace bl = boost::lambda;
          if( _string.compare("partition") != 0 ) return;
          typedef typename T_EVALUATOR :: t_GATraits :: t_Individual t_Individual;
          typedef typename T_EVALUATOR :: t_GATraits :: t_Population t_Population;
          typedef Taboo::Container<t_Individual> t_Taboo;
          typedef boost::function<bool(t_Individual&)> t_Random;
          typedef void( *t_mop )( t_Individual&, size_t, size_t );
         
          _func = bl::bind
                  (
                    &LaDa::GA::Operators::partition_populate< t_Random, t_mop, t_Taboo, t_Population >,
                    bl::protect( bl::bind( &T_EVALUATOR::initialize, bl::var(_evaluator), bl::_1 ) ),
                    &Operators::call_object_mask<t_Individual>,
                    bl::var(_taboo), bl::_1, bl::_2, 50 * bl::_2
                  );
        }
    }
  }// namespace GA
} // namespace LaDa


#endif

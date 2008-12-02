//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_TRUE_CENSUS_H_
#define _LADA_GA_CHECKPOINTS_TRUE_CENSUS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/bind.hpp>

#include <print/stdout.h>
#include <print/xmg.h>
#include <opt/tinyxml.h>

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Prints the number of truly different individuals in the population.
      template< class T_POPULATION >
        void print_truecensus( bool, const T_POPULATION& _pop )
        {
          typename T_POPULATION :: const_iterator i_begin = _pop.begin();
          typename T_POPULATION :: const_iterator i_indiv1 = i_begin;
          typename T_POPULATION :: const_iterator i_end = _pop.end();
          typename T_POPULATION :: const_iterator i_indiv2;
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

      namespace Factory
      {
        //! Factory function for printing the number of truly different individuals in the population.
        template< class T_CHECKPOINT, class T_POPULATION >
          void print_truecensus( T_CHECKPOINT& _checkpoint )
          {
            _checkpoint.connect_statistic( &GA::CheckPoint::print_truecensus<T_POPULATION> ); 
            Print::out << "Will print census at each generation.\n";
            Print::xmg << Print::Xmg::comment << "Will print census at each generation." 
                       << Print::endl;
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::print_truecensus()
        template< class T_CHECKPOINT, class T_POPULATION >
          void (*print_truecensus( const T_CHECKPOINT&, const T_POPULATION& ))
              ( T_CHECKPOINT& ) 
            { return &Factory::print_truecensus< T_CHECKPOINT, T_POPULATION >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

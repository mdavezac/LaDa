//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_AVERAGE_FITNESS_H_
#define _LADA_GA_CHECKPOINTS_AVERAGE_FITNESS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <boost/bind.hpp>

#include <print/stdout.h>
#include <opt/modifiers.h>
#include <opt/tinyxml.h>
#include <opt/types.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Prints the average (scalar) fitness.
      template< class T_POPULATION>
        void print_averagefitness( bool, const T_POPULATION& _pop )
        {
          typename T_POPULATION :: const_iterator i_indiv = _pop.begin();
          typename T_POPULATION :: const_iterator i_indiv_end = _pop.end();
          
          typedef typename T_POPULATION :: value_type 
                                        :: t_IndivTraits 
                                        :: t_ScalarFitness  t_ScalarFitness;
          types::t_real average(0);
          for( ; i_indiv != i_indiv_end; ++i_indiv )
            average += (t_ScalarFitness) i_indiv->fitness();
          
          average /= ( types::t_real ) _pop.size();
          Print::out << "Average Fitness: " << average << "\n\n";
        }

      namespace Factory
      {
        //! Factory function for printing the average (scalar) fitness.
        template< class T_CHECKPOINT, class T_POPULATION >
          void print_averagefitness( T_CHECKPOINT& _checkpoint )
          {
#           ifdef __PGI 
              void (*func)(bool, T_POPULATION const& )
                 = &GA::CheckPoint::print_averagefitness<T_POPULATION>;
              _checkpoint.connect_statistic( func ); 
#           else
              _checkpoint.connect_statistic( &GA::CheckPoint::print_averagefitness<T_POPULATION> ); 
#           endif
            Print::out << "Will print average scalar fitness at each generation.\n";
            Print::xmg << Print::Xmg::comment << "Will print average scalar fitness at each generation." 
                       << Print::endl;
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::print_averagefitness()
        template< class T_CHECKPOINT, class T_POPULATION >
          void (*print_averagefitness( const T_CHECKPOINT&, const T_POPULATION& ))
                  ( T_CHECKPOINT& )
            { return &Factory::print_averagefitness< T_CHECKPOINT, T_POPULATION >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_MAXEVAL_H_
#define _LADA_GA_CHECKPOINTS_MAXEVAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <print/stdout.h>
#include <print/xmg.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Returns true if the T_EVALUATION :: nb_eval < \a _max.
      template< class T_DARWIN >
      bool maxevaluations( const T_DARWIN& _darwin, const types::t_unsigned& _max )
      { 
        if( _darwin.evaluation->nb_eval <= _max ) return true;
        Print::out << "Stopping after " << _darwin.evaluation->nb_eval << " evaluations.\n";
        Print::xmg << Print::Xmg::comment << "Stopping after "
                   << _darwin.evaluation->nb_eval << " evaluations." << Print::endl;
        return false;
      }

      namespace Factory
      {
        //! Factory function for stopping after n evaluations.
        template< class T_CHECKPOINT, class T_DARWIN >
          void maxevaluations( T_CHECKPOINT& _checkpoint,
                               const std::string& _maxeval,
                               const T_DARWIN& _darwin )
          {
            types::t_unsigned maxeval;
            __TRYBEGIN
             maxeval = boost::lexical_cast<types::t_unsigned>( _maxeval );
            __TRYEND(, "Could not parse \"" << _maxeval << "\" into a natural integer.\n" )
            if( maxeval == 0 )
            {
              Print :: xmg << Print::Xmg::comment
                           << "Will not stop on number of evaluations."
                           << Print::endl;
              Print :: out << "Will not stop on number of evaluations.\n";
              return;
            }
            _checkpoint.connect_continuator( boost::bind( &GA::CheckPoint::maxevaluations<T_DARWIN>, 
                                                          boost::cref(_darwin), maxeval ) );
            Print :: xmg << Print :: Xmg :: comment
                         << "Maximum number of evaluations: " << maxeval << "." << Print::endl;
            Print :: out << "Maximum number of evaluations: " << maxeval << ".\n";
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::maxevaluations()
        template< class T_CHECKPOINT, class T_DARWIN >
          void (*maxevaluations( const T_CHECKPOINT&, const T_DARWIN& ))
                  ( T_CHECKPOINT&, const std::string&, const T_DARWIN& )
            { return &Factory::maxevaluations< T_CHECKPOINT, T_DARWIN >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

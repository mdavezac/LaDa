//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_MAXGENERATION_H_
#define _LADA_GA_CHECKPOINTS_MAXGENERATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <print/stdout.h>
#include <print/xmg.h>
#include <opt/debug.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Returns true if the \a _age() <= \a _max.
      inline bool maxgenerations( const GenCount _age, size_t _max )
      {
        if( _age() <= _max ) return true;
        Print::out << "Stopping after " << _age() << " generations.\n";
        Print::xmg << Print::Xmg::comment << "Stopping after "
                   << _age() << " generations." << Print::endl;
        return false;
      }

      namespace Factory
      {
        //! Factory function for stopping after n generations.
        template< class T_CHECKPOINT >
          void maxgenerations( T_CHECKPOINT& _checkpoint,
                               const std::string& _maxgen,
                               const GenCount _age )
          {
            types::t_unsigned maxgen;
            __TRYBEGIN
              maxgen = boost::lexical_cast<types::t_unsigned>(_maxgen);
            __TRYEND(, "Could not parse \"" << _maxgen << "\" into a natural integer.\n" )

            if( maxgen == 0 )
            {
              Print :: xmg << Print::Xmg::comment
                           << "Will not stop on number of generations."
                           << Print::endl;
              Print :: out << "Will not stop on number of generations.\n";
              return;
            }
            _checkpoint.connect_continuator( boost::bind( &GA::CheckPoint::maxgenerations,
                                                          _age, maxgen ) );
            Print :: xmg << Print :: Xmg :: comment
                         << "Maximum number of generations: " << _maxgen << "." << Print::endl;
            Print :: out << "Maximum number of generations: " << _maxgen << ".\n";
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::maxgenerations()
        template< class T_CHECKPOINT >
          void (*maxgenerations( const T_CHECKPOINT& ))
                ( T_CHECKPOINT&, const std::string&, const GenCount )
            { return &Factory::maxgenerations< T_CHECKPOINT >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

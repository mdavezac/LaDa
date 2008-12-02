//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_PRINT_ITERATION_H_
#define _LADA_GA_CHECKPOINTS_PRINT_ITERATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <boost/bind.hpp>

#include <print/xmg.h>
#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Prints iteration #, new results, # of calls to functional.
      template< class T_DARWIN >
        void print_at_iteration( bool _lastcall, bool _eachcall, const T_DARWIN& _darwin )
        {
          if( _lastcall ) { print_lastcall( *_darwin.store, _darwin.counter ); return; }
          if ( not ( _eachcall or _darwin.store->newresults() ) )
          {
            Print::xmg << Print::Xmg::clearall;
            return;
          }
          
          std::string special = "";
          if ( not _darwin.store->newresults() ) special = " ? ";
         
          Print::xmg << Print::Xmg::comment << special << "Iteration " << _darwin.counter() 
                     << Print::endl 
                     << Print::Xmg::comment << special << "Evaluation Calls: " 
                     << _darwin.evaluation->nb_eval << " " << _darwin.evaluation->nb_grad 
                     << Print::endl;
          
          if( _darwin.store->newresults() ) _darwin.store->print_results( _darwin.counter() );
          Print::xmg << Print::flush;
        }

      //! On last call, prints results.
      template< class T_STORE >
        void print_lastcall( const T_STORE& _store, const GenCount& _age )
        {
          Print::xmg << Print::Xmg::comment << "Last Found Result" << Print::endl;
          _store.print_results(_age(), true);
          Print::xmg << Print::flush;
        }

      namespace Factory
      {
        //! Factory function for end-of-generation printing.
        template< class T_CHECKPOINT, class T_DARWIN >
          void print_at_iteration( T_CHECKPOINT& _checkpoints, bool _eachcall,
                                   const T_DARWIN& _darwin )
          {
//           _checkpoints.connect_updater
//           (
//             boost::bind( &GA::CheckPoint::print_at_iteration<T_DARWIN>,
//                          _1, _eachcall, boost::cref(_darwin) )
//           );
          }
      }
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

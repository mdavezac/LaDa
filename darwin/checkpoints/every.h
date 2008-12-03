//
//  Version: $Id$
//
#ifndef _LADA_GA_CHECKPOINTS_EVERY_H_
#define _LADA_GA_CHECKPOINTS_EVERY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/bind.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <tinyxml/tinyxml.h>

#include <print/stdout.h>
#include <print/xmg.h>
#include <mpi/mpi_object.h>
#include <opt/debug.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Calls an updater every \a _every generations.
      template< class T_FUNCTOR >
        void every( bool _lastcall, const GenCount _age, 
                    size_t _every, const T_FUNCTOR _functor )
        {
          if( _every == 0 ) { _functor(); return; }
          if( not ( _age() == 0 or _age() % _every == 0 or _lastcall ) ) return;
          _functor();
        }

      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::max_generations()
        template< class T_FUNCTOR >
          void (*every( const T_FUNCTOR& ))
                  ( bool, const GenCount, size_t, const T_FUNCTOR )
            { return &Factory::every< T_FUNCTOR >; }
      }
  

    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_STOP_ONFILE_H_
#define _LADA_GA_CHECKPOINTS_STOP_ONFILE_H_

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
      //! Returns true if \a _path does not exist.
      inline bool stop_onfile( const boost::filesystem::path _path )
      {
        __ROOTCODE
        (
          (*LaDa::mpi::main),
          if( not boost::filesystem::exists( _path ) ) return true;
          Print::out << "Stopping on finding file " << _path << ".\n";
          Print::xmg << Print::Xmg::comment << "Stopping on finding file "
                     << _path << "." << Print::endl;
          return false;
          
        )
        __MPICODE( return true; )
      }

      namespace Factory
      {
        //! Factory function for stopping on finding file.
        template< class T_CHECKPOINT >
          void stop_onfile( T_CHECKPOINT& _checkpoint, const std::string& _path )
          {
            const boost::filesystem::path path( _path );
            __ROOTCODE
            (
              (*LaDa::mpi::main),
              if( boost::filesystem::exists( path ) )
              __BEGINGROUP__
                __DOASSERT( not boost::filesystem::is_regular( path ),
                            "Stop file " << path << " with unknown (eg not regular file)"
                            " format found on startup.\n" )
                boost::filesystem::remove_all( path );
                std::cerr << "Found and deleted stop file " << path << " on startup.\n";
              __ENDGROUP__
            )
            _checkpoint.connect_continuator( boost::bind( &GA::CheckPoint::stop_onfile, path ) );
            Print :: xmg << Print :: Xmg :: comment
                         << "Will stop on finding file " << path << "." << Print::endl;
            Print :: out << "Will stop on finding file " << path << ".\n";
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::max_generations()
        template< class T_CHECKPOINT >
          void (*stop_onfile( const T_CHECKPOINT& ))( T_CHECKPOINT&, const std::string& ) 
            { return &Factory::stop_onfile< T_CHECKPOINT >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

//
//  Version: $Id: print_callbacks.h 871 2008-12-02 05:24:15Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <print/stdout.h>
#include <print/xmg.h>
#include <opt/tinyxml.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "init_outputfiles.h"

namespace LaDa
{

  namespace GA
  {
    void init_outputfiles( const boost::filesystem::path &_input )
    {
      __TRYBEGIN
      TiXmlDocument doc;
      TiXmlHandle docHandle( &doc ); 
      opt::read_xmlfile( _input, doc );

      boost::filesystem::path xmgpath( "xmg" );
      boost::filesystem::path outpath( "out" );

      const TiXmlElement *child = docHandle.FirstChild("Job")
                                           .FirstChild("GA" )
                                           .FirstChild("Filenames").Element();
      for(; child; child = child->NextSiblingElement() )
        if( child->Attribute("out") )
          outpath = child->Attribute("out");
        else if( child->Attribute("xmgrace") )
          xmgpath = child->Attribute("xmgrace");

      __MPISERIALCODE( 
        // MPI code
        __TRYCODE(    Print::out.sync_filename( outpath );
                      Print::xmg.sync_filename( xmgpath );,
                   "Caught error while synchronizing output filenames\n" 
        ),
        // Serial code
        Print::xmg.init( xmg_filename );
        Print::out.init( out_filename );
      )
      __TRYEND(,"Could not start output files.\n")
    }
  } // namespace GA
} // namespace LaDa


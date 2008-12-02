//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_XYZANIM_H_
#define _LADA_GA_CHECKPOINTS_XYZANIM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/path.hpp>
#include <fstream>
#include <sstream>
#include <crystal/structure.h>

#include "apply2best.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Prints an individual to a file in XYZ format.
      template< class T_INDIVIDUAL >
        void xyz_anim( const T_INDIVIDUAL& _individual, const std::string& _filename, 
                       const GenCount& _age, const Crystal::Structure& _structure )
        {
          std::ofstream file( _filename.c_str(), std::ios_base::app | std::ios_base::out); 
          if( file.fail() ) 
          {
            std::cerr << "Could not open " << _filename << " for writing.\n"
                      << "Will not print xyz animation.\n";
            return;
          }
          file << "# Iteration: " << _age() << "  --  " << _individual.fitness();
          Crystal :: Structure structure( _structure );
          structure << _individual.Object();
          structure.print_xyz( file );
          file.close();
        };
      
      namespace Factory
      {
        //! A factory for creating an xyz animation of best result at each generation.
        template< class T_CHECKPOINT, class T_DARWIN >
          void xyz_anim( T_CHECKPOINT& _checkpoint,
                              const std::string &_filename,
                              const T_DARWIN& _darwin, 
                              const Crystal::Structure& _structure )
          {
            typedef typename T_DARWIN :: t_GATraits :: t_Individual t_Individual;
            boost::filesystem::path path( _filename );
            _checkpoint.connect_updater
            (
              apply2best
              ( 
                 _darwin, 
                 boost::bind
                 ( 
                   &GA::CheckPoint::xyz_anim< t_Individual >,
                   _1,  path.string(), boost::cref( _darwin.counter ), boost::cref( _structure ) 
                 ) 
              )
            );
          }
      } // namespace factory.
      //! Collects helper functions for obtaining function address for binding in lambdas.
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::xyz_anim()
        template< class T_CHECKPOINT, class T_DARWIN >
          void (*xyz_anim( const T_CHECKPOINT&, const T_DARWIN& ))
            ( T_CHECKPOINT&, const std::string&, const T_DARWIN&, const Crystal::Structure& )
              { return &Factory::xyz_anim< T_CHECKPOINT, T_DARWIN >; }
      } // namespace AddressOf
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

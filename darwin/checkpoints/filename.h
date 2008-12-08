//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_FILENAME_H_
#define _LADA_GA_CHECKPOINTS_FILENAME_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <boost/filesystem/path.hpp>

#include <print/stdout.h>


namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      namespace Factory
      {
        //! Factory function for stopping after n evaluations.
        template< class T_CHECKPOINT >
          void filename( T_CHECKPOINT&,
                         const std::string& _att,
                         const std::string& _name,
                         boost::filesystem::path& _result )
          {
            Print::out << _name << _att << "\n";
            _result = _att;
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::maxevaluations()
        template< class T_CHECKPOINT >
          void (*filename( const T_CHECKPOINT& ))
                  ( T_CHECKPOINT&, const std::string&, const std::string&, boost::filesystem::path& )
            { return &Factory::filename< T_CHECKPOINT >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif

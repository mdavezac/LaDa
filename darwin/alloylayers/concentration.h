//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_CONCENTRATION_H_
#define _DARWIN_ALLOY_LAYERS_CONCENTRATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/serialization/string.hpp>
#include <boost/lexical_cast.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>

namespace LaDa
{
  namespace GA
  {

    namespace Keepers
    {
      //! Keeps track of concentration.
      class Concentration
      {
        friend class boost::serialization::access;
        public:
          //! Concentration on each site 0.
          types::t_real x;
          //! Concentration on each site 1.
          types::t_real y;
          //! Constructor.
          Concentration() : x(-2), y(-2) {}
          //! Copy Constructor.
          Concentration( const Concentration &_c ) : x( _c.x ), y( _c.y ) {}
          //! Loads concentration from XML.
          bool Load ( const TiXmlElement &_node );
          //! Saves concentration to XML.
          bool Save ( TiXmlElement &_node ) const;
        private:
          //! Serializes concentration.
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
            { _ar & x; _ar & y; }
      };

    } // end of Keepers namespace.

  } // namespace GA
} // namespace LaDa
#endif

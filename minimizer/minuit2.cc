#include "LaDaConfig.h"

#include <iomanip>
#include <boost/lexical_cast.hpp>

#include "minuit2.h"

namespace LaDa
{
  namespace Minimizer {

    LADA_REGISTER_MINIMIZER_VARIANT_SOURCE( Minuit2, "Minuit2 Migrad" )


    bool Minuit2 :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = opt::find_node( _node, "Minimizer", "type", "minuit2" );
      if( not parent ) return false;
      if( parent->Attribute( "maxcalls" ) )
        itermax = boost::lexical_cast< types::t_unsigned >( parent->Attribute("maxcalls") );
      else if( parent->Attribute( "itermax" ) )
        itermax = boost::lexical_cast< types::t_unsigned >( parent->Attribute("itermax") );
      if( parent->Attribute("tolerance") ) 
        tolerance = boost::lexical_cast< types::t_real >( parent->Attribute("tolerance") );
      if( parent->Attribute("strategy") ) 
      {
        const std::string value = parent->Attribute("strategy");
        if( value == "fast" ) strategy = 0u;
        else if( value == "slow" ) strategy = 1u;
        else if( value == "slowest" ) strategy = 2u;
        else __THROW_ERROR( "Unknown strategy " + value + "\n" );
      }
      if( parent->Attribute("uncertainties") ) 
        uncertainties = boost::lexical_cast< types::t_real >( parent->Attribute("uncertainties") );
      if( parent->Attribute("verbose") ) 
      {
        const std::string value( parent->Attribute("verbose") );
        if( value == "true" or value == "TRUE" or value == "T" or value == "t" ) 
          verbose = true;
        else if( value == "false" or value == "FALSE" or value == "F" or value == "f" ) 
          verbose = false;
        else verbose = boost::lexical_cast<bool>( parent->Attribute("verbose") );
      }
      if( parent->Attribute("up") ) 
        up = boost::lexical_cast< types::t_real >( parent->Attribute("uncertainties") );
      if( parent->Attribute("gradient") ) 
      {
        const std::string value( parent->Attribute("gradient") );
        if( value == "true" or value == "TRUE" or value == "T" or value == "t" ) 
          use_gradient = true;
        else if( value == "false" or value == "FALSE" or value == "F" or value == "f" ) 
          use_gradient = false;
        else verbose = boost::lexical_cast<bool>( parent->Attribute("gradient") );
      }
      return true;
    }
  }
} // namespace LaDa

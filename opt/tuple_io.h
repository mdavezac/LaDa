#ifndef _OPT_TUPLES_H_
#define _OPT_TUPLES_H_

#include "LaDaConfig.h"

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include "types.h"
#include "debug.h"

namespace LaDa
{
  namespace opt
  {
    namespace tuples
    {
      template< class T_TYPE >
        bool read( const std::string &_string, 
                   boost::tuple<T_TYPE, T_TYPE, T_TYPE>& _t )
        {
          namespace bt = boost::tuples;
          __TRYBEGIN
            boost::regex re("\\(\\s*(\\S+)\\s*,?"
                            "\\s*(\\S+)\\s*,?"
                            "\\s*(\\S+)\\s*\\)");
            boost::match_results<std::string::const_iterator> what;
            if( not boost::regex_search( _string, what, re ) )
            {
              re = boost::regex("(\\S+)\\s*,?"
                                "\\s*(\\S+)\\s*,?"
                                "\\s*(\\S+)");
              if( not boost::regex_search( _string, what, re ) )  return false;
            }
            bt::get<0>(_t) = boost::lexical_cast<T_TYPE>( what.str(1) );
            bt::get<1>(_t) = boost::lexical_cast<T_TYPE>( what.str(2) );
            bt::get<2>(_t) = boost::lexical_cast<T_TYPE>( what.str(3) );
          __TRYEND(,"Could not parse tuple. Expects something like \"( a b c )\".\n")
          return true;
        }
    }
  }
} // namespace LaDa

#endif

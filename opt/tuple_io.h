//
//  Version: $Id$
//
#ifndef _OPT_TUPLES_H_
#define _OPT_TUPLES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include "types.h"

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
          const boost::regex re("\\(?(\\s+)?(\\S+)(\\s+|,)(\\S+)(\\s+|,)(\\S+)(\\s+)?\\)?");
          boost::match_results<std::string::const_iterator> what;
          if( not boost::regex_search( _string, what, re ) ) return false;
          bt::get<0>(_t) = boost::lexical_cast<types::t_int>( what.str(2) );
          bt::get<1>(_t) = boost::lexical_cast<types::t_int>( what.str(4) );
          bt::get<2>(_t) = boost::lexical_cast<types::t_int>( what.str(6) );
        __TRYEND(,"Could not parse tuple.\n")
        return true;
      }
  }
}

#endif

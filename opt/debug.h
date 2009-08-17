//
//  Version: $Id$
//

#ifndef _OPT_DEBUG_H_
#define _OPT_DEBUG_H_
#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <iostream>

#if defined(foreach)
#  warning "foreach macro already exists"
#endif
#define foreach BOOST_FOREACH

#define __SPOT_ERROR __FILE__ << ", line: " << __LINE__ << "\n" 
#define __CATCHCODE(code, error)\
        catch(...)\
        {\
          code;\
          std::cerr << __SPOT_ERROR << error; \
          throw 0; \
        }
#define __THROW_ERROR(error) \
        {\
          std::cerr << __SPOT_ERROR << error ;\
          throw 0; \
        }
#define __TRYCODE(code,error) try { code } \
        catch ( std::exception &_e )\
        __THROW_ERROR( error << _e.what() )
#define __ASSERTCATCHCODE(condition, code, error) \
          if( condition ) { code; __THROW_ERROR( error ) }
#define __DOASSERT( condition, error ) \
          __ASSERTCATCHCODE( condition, ,error )
#define __TRYASSERT(condition, error) \
          __TRYCODE( __DOASSERT( condition, "" ), error )
#define LADA_ASSERT(a,b) __ASSERT( (not (a) ), b)
#define __TRYBEGIN try {
#define __TRYEND( code, error ) } __CATCHCODE( code, error ) 
#define __COMMA__ ,
#define __BEGINGROUP__ {
#define __ENDGROUP__ }


#ifdef _LADADEBUG
# define LADA_DEBUG
# define __DODEBUGCODE(code) code
# define __TRYDEBUGCODE(code, error) __TRYCODE(code, error)
# define __DOTRYDEBUGCODE(code, error) __TRYCODE(code, error)
# define __ASSERT( condition, error ) __DOASSERT(condition, error)
# define __NDEBUGCODE( code ) 
# define __DEBUGCATCHCODE( code, error ) __CATCHCODE( code, error )
# define __DEBUGTRYBEGIN try {
# define __DEBUGTRYEND( code, error ) } __CATCHCODE( code, error ) 
#else
# define __ASSERT( condition, error ) 
# define __DODEBUGCODE(code) 
# define __TRYDEBUGCODE(code, error) code
# define __DOTRYDEBUGCODE(code, error) 
# define __NDEBUGCODE( code ) code
# define __DEBUGCATCHCODE( code, error ) 
# define __DEBUGTRYBEGIN 
# define __DEBUGTRYEND( code, error ) 
#endif

#endif

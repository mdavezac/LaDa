//
//  Version: $Id$
//

#ifndef _OPT_DEBUG_H_
#define _OPT_DEBUG_H_
#include <iostream>
#include <stdexcept>
#include <boost/foreach.hpp>

#if defined(foreach)
#  warning "foreach macro already exists"
#endif
#define foreach BOOST_FOREACH

#define __SPOT_ERROR __FILE__ << ", line: " << __LINE__ << "\n" 
#define __CATCHCODE(code, error)\
        catch( std::exception &_e )\
        {\
          code;\
          std::ostringstream sstr;\
          sstr << __SPOT_ERROR << error << _e.what();\
          throw std::runtime_error( sstr.str() );\
        }
#define __THROW_ERROR(error) \
        {\
          std::ostringstream sstr;\
          sstr << __SPOT_ERROR << error ;\
          throw std::runtime_error( sstr.str() );\
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
#define __TRYBEGIN try {
#define __TRYEND( code, error ) } __CATCHCODE( code, error ) 

#define __CATCHBPO \
  } \
  catch ( boost::program_options::invalid_command_line_syntax &_b)\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Something wrong with the command-line input.\n"\
              << _b.what() << std::endl;\
    return 0;\
  }\
  catch ( boost::program_options::invalid_option_value &_i )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Argument of option in command-line is invalid.\n"\
              << _i.what() << std::endl;\
    return 0;\
  }\
  catch ( boost::program_options::unknown_option &_u)\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Unknown option in command-line.\n"\
              << _u.what() << std::endl;\
    return 0;\
  }\
  catch (  boost::program_options::too_many_positional_options_error &_e )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"\
              << "Too many arguments in command-line.\n"\
              << _e.what() << std::endl;\
    return 0;\
  }\
  catch ( std::exception &e )\
  {\
    std::cout << "Caught error while running " << __PROGNAME__ \
              << "\n" << e.what() << std::endl;\
    return 0;\
  }

#ifdef _LADADEBUG
#define __DODEBUGCODE(code) code
#define __TRYDEBUGCODE(code, error) __TRYCODE(code, error)
#define __DOTRYDEBUGCODE(code, error) __TRYCODE(code, error)
#define __ASSERT( condition, error ) __DOASSERT(condition, error)
#define __NDEBUGCODE( code ) 
#define __DEBUGCATCHCODE( code, error ) __CATCHCODE( code, error )
#define __DEBUGTRYBEGIN try {
#define __DEBUGTRYEND( code, error ) } __CATCHCODE( code, error ) 
#else
#define __ASSERT( condition, error ) 
#define __DODEBUGCODE(code) 
#define __TRYDEBUGCODE(code, error) code
#define __DOTRYDEBUGCODE(code, error) 
#define __NDEBUGCODE( code ) code
#define __DEBUGCATCHCODE( code, error ) 
#define __DEBUGTRYBEGIN 
#define __DEBUGTRYEND( code, error ) 
#endif

#endif

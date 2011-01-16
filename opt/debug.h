#ifndef _OPT_DEBUG_H_
#define _OPT_DEBUG_H_
#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <iostream>
#include <boost/throw_exception.hpp>

#if defined(foreach)
#  warning "foreach macro already exists"
#endif
#define foreach BOOST_FOREACH

#define LADA_SPOT_ERROR __FILE__ << ", line: " << __LINE__ << "\n" 
#define LADA_CATCHCODE(code, error)\
        catch(...)\
        {\
          code;\
          std::cerr << LADA_SPOT_ERROR << error; \
          throw 0; \
        }
#define LADA_THROW_ERROR(error) \
        {\
          std::cerr << LADA_SPOT_ERROR << error ;\
          throw 0; \
        }
#define LADA_TRY_CODE(code,error) try { code } \
        catch ( std::exception &_e )\
        LADA_THROW_ERROR( error << _e.what() )
#define LADA_NASSERT_CATCHCODE(condition, code, error) \
          if( condition ) { code; LADA_THROW_ERROR( error ) }
#define LADA_DO_NASSERT( condition, error ) \
          LADA_NASSERT_CATCHCODE( condition, ,error )
#define LADA_TRY_ASSERT(condition, error) \
          LADA_TRY_CODE( LADA_DO_NASSERT( condition, "" ), error )
#define LADA_ASSERT(a,b) LADA_NASSERT( (not (a) ), b)
#define LADA_DOASSERT(a,b) LADA_DO_NASSERT( (not (a) ), b)
#define LADA_TRY_BEGIN try {
#define LADA_TRY_END( code, error ) } LADA_CATCHCODE( code, error ) 
#define LADA_COMMA ,
#define LADA_BEGINGROUP {
#define LADA_ENDGROUP }

#ifdef LADA_DEBUG
# include <boost/throw_exception.hpp>
# define LADA_BASSERT(condition, error) if(not condition) BOOST_THROW_EXCEPTION(error)
# define LADA_NASSERT(condition, error ) LADA_DO_NASSERT(condition, error)
# define LADA_DEBUG_TRY_BEGIN try {
# define LADA_DEBUG_TRY_END( code, error ) } LADA_CATCHCODE( code, error ) 
#else
# define LADA_BASSERT(condition, error) 
# define LADA_NASSERT( condition, error ) 
# define LADA_DEBUG_TRY_BEGIN 
# define LADA_DEBUG_TRY_END( code, error ) 
#endif

#endif

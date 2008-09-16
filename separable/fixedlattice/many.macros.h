//
//  Version: $Id$
//
#ifndef _SEPARABLES_FIXEDLATTICES_MANY_MACROS_H_
#define _SEPARABLES_FIXEDLATTICES_MANY_MACROS_H_

# ifdef VOID_MFUNC
#  error VOID_MFUNC already defined.
# endif 
# ifdef VOID_MFUNC_DECLARE
#  error VOID_MFUNC_DECLARE already defined.
# endif 
# define VOID_MFUNC( FunctionName ) Void_Unary_Apply ## FunctionName
# define VOID_MFUNC_DECLARE( FunctionName ) \
   struct VOID_MFUNC( FunctionName )\
   {\
     typedef void result_type;\
     template< class T > result_type operator()( T &_t ) \
     {\
       foreach( const typename T::value_type &_val, _t )\
         _val.FunctionName(); \
     }\
   }; 

# ifdef ACC_MFUNC
#  error ACC_MFUNC already defined.
# endif 
# ifdef ACC_MFUNC_DECLARE
#  error ACC_MFUNC_DECLARE already defined.
# endif 
# define ACC_MFUNC( FunctionName ) AccApply ## FunctionName
# define ACC_MFUNC_DECLARE( FunctionName, Type ) \
   struct ACC_MFUNC( FunctionName )\
   {\
     typedef Type result_type;\
     template< class T > result_type operator()( T &_t, result_type _r )\
     {\
       foreach( const typename T::value_type &_val, _t )\
         _r += _val.FunctionName(); \
       return _r;\
     }\
   }; 

# ifdef U_MFUNC
#  error U_MFUNC already defined.
# endif 
# ifdef U_MFUNC_DECLARE
#  error U_MFUNC_DECLARE already defined.
# endif 
# define U_MFUNC( FunctionName ) UApply ## FunctionName
# define U_MFUNC_DECLARE(FunctionName,Code,ArgType) \
   struct U_MFUNC(FunctionName) \
   { \
     typedef void result_type;\
     ArgType arg;\
     U_MFUNC( FunctionName ) ( const ArgType _t ) : arg( _t ) {}\
     U_MFUNC( FunctionName ) ( const U_MFUNC( FunctionName ) &_c ) : arg( _c.arg ) {}\
     template< class T > result_type operator()( T &_t )\
     {\
       foreach( const typename T::value_type &_val, _t )\
         _val Code;\
     }\
   };

# if defined( MANYHEAD ) || defined(INMANY) || defined(INMANY2)
#   error "Macros with same names."
# endif
# define MANYHEAD \
    Many<T_TRAITS> 
#  define INMANY( var ) \
     template< class T_TRAITS >  var MANYHEAD
#  define INMANY2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 MANYHEAD


#else

# ifdef VOID_MFUNC
#   undef VOID_MFUNC 
# endif 
# ifdef VOID_MFUNC_DECLARE
#   undef VOID_MFUNC_DECLARE
# endif 
# ifdef ACC_MFUNC
#   undef ACC_MFUNC 
# endif 
# ifdef ACC_MFUNC_DECLARE
#   undef ACC_MFUNC_DECLARE
# endif 
# ifdef U_MFUNC
#   undef U_MFUNC 
# endif 
# ifdef U_MFUNC_DECLARE
#   undef U_MFUNC_DECLARE
# endif 
# ifdef MANYHEAD
#   undef MANYHEAD
# endif
# ifdef INMANY
#  undef INMANY
# endif
# ifdef INMANY2
#  undef INMANY2
# endif


#endif

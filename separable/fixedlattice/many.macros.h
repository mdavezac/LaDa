//
//  Version: $Id$
//
#ifndef _SEPARABLES_FIXEDLATTICES_MANY_MACROS_H_
#define _SEPARABLES_FIXEDLATTICES_MANY_MACROS_H_


#ifdef PHOENIX_MEMFUNC
#  error PHOENIX_MEMFUNC already defined.
# endif 
# ifdef PHOENIX_MEMFUNC_DECLARE
#  error PHOENIX_MEMFUNC_DECLARE already defined.
# endif 
#define PHOENIX_MEMFUNC( memfunc ) memfunc_ ## memfunc ## _impl
#define PHOENIX_MEMFUNC_DECLARE( memfunc, resulttype ) \
  struct PHOENIX_MEMFUNC( memfunc ) \
  { \
      template <typename Arg> \
      struct result \
      {  \
        typedef resulttype type; \
      }; \
      \
      template <typename Arg> \
      typename boost::disable_if\
               < boost::is_same< typename result<Arg> :: type, void >, \
                 typename result<Arg> :: type > :: type \
        operator()(Arg& n) const  {  return n.memfunc();  } \
      template <typename Arg> \
      typename boost::enable_if \
               < boost::is_same< typename result<Arg> :: type, void >, void > :: type\
        operator()(Arg& n) const  {  n.memfunc();  } \
  }; \

#ifdef PHOENIX_MEMFUNC1
#  error PHOENIX_MEMFUNC1 already defined.
# endif 
# ifdef PHOENIX_MEMFUNC1_DECLARE
#  error PHOENIX_MEMFUNC1_DECLARE already defined.
# endif 
#define PHOENIX_MEMFUNC1( memfunc ) memfunc1_ ## memfunc ## _impl
#define PHOENIX_MEMFUNC1_DECLARE( memfunc, resulttype ) \
  struct PHOENIX_MEMFUNC1( memfunc ) \
  { \
      template <typename Arg, typename Arg1 = void > \
      struct result \
      {  \
        typedef resulttype type; \
      }; \
      \
      template <typename Arg, typename Arg1> \
      typename boost::disable_if\
               < boost::is_same< typename result<Arg, Arg1> :: type, void >, \
                 typename result<Arg, Arg1> :: type > :: type \
        operator()(Arg& n, Arg1& a1 ) const  {  return n.memfunc( a1 );  } \
      template <typename Arg, typename Arg1 > \
      typename boost::enable_if \
               < boost::is_same< typename result<Arg, Arg1> :: type, void >, void > :: type\
        operator()(Arg& n, Arg1& a1) const  {  n.memfunc( a1 );  } \
  }; \

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
       foreach( typename T::const_reference &_val, _t )\
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
     template< class T > result_type operator()( T &_t, result_type result )\
     {\
       foreach( typename T::const_reference &_val, _t )\
         result += _val.FunctionName(); \
       return result;\
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
     U_MFUNC( FunctionName ) ( ArgType _t ) : arg( _t ) {}\
     U_MFUNC( FunctionName ) ( const U_MFUNC( FunctionName ) &_c ) : arg( _c.arg ) {}\
     template< class T > result_type operator()( T &_t )\
     {\
       foreach( typename T::const_reference &_val, _t )\
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
# ifdef PHOENIX_MEMFUNC
#  undef PHOENIX_MEMFUNC
# endif
# ifdef PHOENIX_MEMFUNC_DECLARE
#  undef PHOENIX_MEMFUNC_DECLARE
# ifdef PHOENIX_MEMFUNC1
#  undef PHOENIX_MEMFUNC1
# endif
# ifdef PHOENIX_MEMFUNC_DECLARE1
#  undef PHOENIX_MEMFUNC_DECLARE1
# endif
# endif


#endif

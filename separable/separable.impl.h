//
//  Version: $Id$
//

#include <boost/xpressive/algorithms.hpp>
#include <boost/xpressive/primitive.hpp>
#include <opt/algorithms.h>

namespace LaDa
{
  namespace Separable
  {

    template< class T_B, template<class> class T_G, template<class> class T_S >
      template< class T_ITERATOR >
        typename Base<T_B, T_G, T_S>::t_Return
          Base<T_B,T_G,T_S> :: operator()( T_ITERATOR _first,
                                           T_ITERATOR _last ) const
          {
            typedef typename t_Basis :: const_iterator function_iterator;
            typedef typename t_Coefs :: const_iterator coef_iterator;
            __ASSERT( basis.size() != coefs.size(),
                      "Incoherent basis/coefficient sizes.\n" )
            __ASSERT( basis.size() != ( _last - _first ),
                      "Incoherent basis/coefficient sizes.\n" )
            function_iterator i_func = basis.begin();
            coef_iterator i_coef = coefs.begin();
            t_Return result( scalarop( *i_coef, (*i_func)( *_first ) ) );
            for( ++_first, ++i_func, ++i_coef; 
                 _first != _last; ++i_func, ++_first, ++i_coef ) 
              result = groupop( result, scalarop( *i_coef, (*i_func)( *_first ) ) ); 
            return result; 
          }

    template< class T_B, template<class> class T_G, template<class> class T_S >
      typename Base<T_B,T_G,T_S>::t_Return
        Base<T_B,T_G,T_S> :: operator()( t_Arg _arg ) const
        {
          typedef typename t_Basis :: const_iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT( basis.size() != coefs.size(),
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          function_iterator i_func_end = basis.end();
          coef_iterator i_coef = coefs.begin();
          t_Return result( scalarop( *i_coef, (*i_func)( _arg ) ) );
          for( ++i_func, ++i_coef; i_func != i_func_end; ++i_func, ++i_coef )
            result = groupop( result, scalarop( *i_coef, (*i_func)( _arg ) ) );
          return result;
        }

    template< class T_B, template<class> class T_G, template<class> class T_S >
      template< class T_ITERATOR >
        void Base<T_B,T_G,T_S> :: set( T_ITERATOR _first, T_ITERATOR _last )
        {
          __ASSERT( (_last-_first) != coefs.size(), 
                    "Too many values in range.\n" )
          std::copy( _first, _last, coefs.begin() );
        }

    template< class T_B, template<class> class T_G, template<class> class T_S >
      template<class ARCHIVE>
        void Base<T_B, T_G, T_S> :: serialize( ARCHIVE & _ar,
                                               const unsigned int _version)
        {
          _ar & basis;
          _ar & coefs;
        }

    template< class T_BASIS > template< class T_ITERATOR >
      typename Function<T_BASIS>::t_Return
        Function<T_BASIS> :: operator()( T_ITERATOR _first,
                                         T_ITERATOR _last ) const
        {
          typedef typename t_Basis :: const_iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() != coefs.size(),
                      "Incoherent basis/coefficient sizes: "
                    << basis.size() << "!=" << coefs.size() << ".\n" )
          function_iterator i_func = basis.begin();
          function_iterator i_func_end = basis.end();
          coef_iterator i_coef = coefs.begin();
          t_Return result( scalarop( *i_coef, (*i_func)( _first, _last ) ) );
          for( ++i_func, ++i_coef; i_func != i_func_end; ++i_func, ++i_coef ) 
            result = groupop( result, 
                              scalarop( *i_coef, (*i_func)( _first, _last ) ) ); 
          return result; 
        }
          
    template< class T_B, template<class> class T_G, template<class> class T_S>
    std::ostream& operator<<( std::ostream& _stream,
                              const Base<T_B, T_G, T_S>& _func )
    {
      namespace bl = boost::lambda;
      if( _func.basis.empty() ) return _stream;

      std::ostringstream sstr;
      opt::concurrent_loop
      (
        _func.basis.begin(), _func.basis.end(), _func.coefs.begin(),
        bl::var( sstr ) << bl::_2 << bl::constant(" ")
                        << bl::_1 << bl::constant("\n") 
      );
      std::string str = sstr.str();
      size_t i = str.rfind('\n');
      if( i != std::string::npos ) str = str.substr( 0, i );
      _stream << _func.name << bx::regex_replace( ":\n" + str, bx::_n, "\n   ");
      return _stream;
    }

    //! \cond
    namespace details
    { 
      //! Wraps an unary free function.
      template< class T_FUNCTION >
         class FreeFunction< T_FUNCTION, 1 >
         {
           //! Traits of the function
           typedef boost::function_traits<T_FUNCTION> t_Traits;
           public:
             //! Type of the free function.
             typedef T_FUNCTION t_Function;
             //! Type of the return.
             typedef typename t_Traits :: result_type t_Return;
             //! Type of the argument.
             typedef typename t_Traits :: arg1_type t_Arg;

             //! Initializes both zero order and gradient function pointers.
             FreeFunction   ( t_Function* _func, t_Function *_grad = NULL) 
                          : func( _func ), grad( _grad ){}
             //! Returns the zero order evaluation.
             t_Return operator()( t_Arg &_x ) { return (*func)(_x); }
             //! Returns the gradient evaluation.
             t_Return gradient( t_Arg &_x ) { return (*grad)(_x); }

           protected:
             //! Function pointer to zero order.
             t_Function *func;
             //! Function pointer to gradient function.
             t_Function *grad;
         };
    } // end namespace details
    //! \endcond


  }
} // namspace LaDa

//
//  Version: $Id$
//

namespace Separable
{
  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITERATOR >
      typename Base<T_BASIS, T_LAW >::t_Return
        Base<T_BASIS, T_LAW> :: operator()( T_ITERATOR _first, T_ITERATOR _last ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() == coefs.size()
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( *_first ) ) ) ;
          for( ++_first, ++i_func, ++i_coef; 
               _first != _last; ++i_func, ++_first, ++i_coef )
            groupop( result, scalarop( *i_coef, (*i_func)( *_first ) ) );
          return result;
        }
  template< class T_BASIS, template<class> class T_LAW >
      typename Base<T_BASIS, T_LAW>::t_Return
        Base<T_BASIS, T_LAW> :: operator()( t_Arg _arg ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT( basis.size() == coefs.size(),
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          function_iterator i_func_end = basis.end();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( _arg ) ) ) ;
          for( ++i_func, ++i_coef; 
               i_func != i_func_end; ++i_func, ++i_coef )
            groupop( result, scalarop( *i_coef, (*i_func)( _arg ) ) );
          return result;
        }
  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITERATOR >
      typename Base<T_BASIS, T_LAW>::t_Return
        Base<T_BASIS, T_LAW> :: set( T_ITERATOR _first, T_ITERATOR _last )
        {
          __ASSERT( (_last-_first) == coefs.size(), 
                    "Too many values in range.\n" )
          std::copy( _first, _last, coefs.begin() );
        }
  template< class T_BASIS, template<class> class T_LAW >
    template<class ARCHIVE>
      void BASE<T_BASIS, T_LAW> :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        ar & basis;
        ar & coefs;
      }

  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITERATOR >
      typename Function<T_BASIS, T_LAW >::t_Return
        Function<T_BASIS, T_LAW> :: operator()( T_ITERATOR _first, T_ITERATOR _last ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() == coefs.size()
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( _first ) ) ) ;
          for( ++i_func, ++i_coef; 
               _first != _last; ++i_func, ++i_coef )
            groupop( result, scalarop( *i_coef, (*i_func)( _first ) ) );
          return result;
        }
  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ARGIT, class T_RETIT >
      typename Function<T_BASIS, T_LAW >::t_Return
        Function<T_BASIS, T_LAW> :: gradient(T_ARGIT _in, T_RETIT _out ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() == coefs.size()
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          function_iterator i_func_end = basis.end();
          coef_iterator i_coef = coefs.begin();
          std::vector< t_Return > store( basis.size() );
          for(; i_fun != i_fun_end; ++i_func, ++i_coef )
          {
            std::fill( store.begin(), store.end(), t_Coef(0) );
            ( i_func->gradient( _in, store.begin() ) );
            std::transform( store.begin(), store.end(), _out, _out
                              boost::lambda::_1 * boost::constant( *i_coef )
                            + boost::lambda::_2 );
          }
        }

  //! \cond
  namespace details
  { 

    template< class T_BASE >
      typename T_RETIT::value_type gradient( T_BASE &_base,
                                             typename T_BASE :: t_Arg _arg )
      {
        return Gradient< typename T_BASE::t_Law >()
                       ( _base.coefs.begin(), _base.coef.end(), 
                         _base.basis.begin(), _arg );
      }
    template< class T_BASE, class T_ARGSIT, class T_RETIT >
      typename T_RETIT::value_type gradient( T_BASE &_base,
                                             T_ARGSIT _arg,
                                             T_RETIT _ret  )
      {
        return Gradient< typename T_BASE::t_Law >()
                       ( _base.coefs.begin(), _base.coef.end(), 
                         _base.basis.begin(), _arg, _ret );
      }
    
    
    template< template<class> class T_LAW > 
    class Gradient
    {
      template< class T_COEF, class T_FUNC, class T_ARG >
        T_COEF operator( T_COEF _first, T_COEF _last, 
                         T_FUNC _func, T_ARG _arg ) 
        { throw std::runtime_error( "Not implemented." ); }
      template< class T_COEF, class T_FUNC, class T_ARG, class T_RET>
        T_COEF operator( T_COEF _first, T_COEF _last, 
                         T_FUNC _func, T_ARG _arg, T_RET _ret ) 
        { throw std::runtime_error( "Not implemented." ); }
    };
    template<> class Gradient< std::plus >
    {
      template< class T_COEF, class T_FUNC, class T_ARG >
        T_COEF operator( T_COEF _first, T_COEF _last, 
                         T_FUNC _func, T_ARG _arg ) 
        {
          T_COEF result(0);
          for( _first ! = _last; ++_coef, ++_func )
            result += (*_coef ) * ( _func->gradient( _arg ) );
          return result;
        }
      template< class T_COEF, class T_FUNC, class T_ARG, class T_RET>
        T_COEF operator( T_COEF _first, T_COEF _last, 
                         T_FUNC _func, T_ARG _arg, T_RET _ret ) 
        {
          for( _first ! = _last; ++_coef, ++_args, ++_func, ++_args, ++_ret)
            *_ret = (*_coef ) * ( _func->gradient(*_args ) );
          return result;
        }
    };
    template<> class Gradient< std::multiplies >
    {
      template< class T_COEF, class T_FUNC, class T_ARG >
        T_COEF operator( T_COEF _first, T_COEF _last, 
                         T_FUNC _func, T_ARG _arg ) 
        {
          // Stores evaluations and gradients.
          typedef std::pair< T_COEF, T_COEF > t_Pair;
          typedef std::vector< t_Pair > t_Store;
          t_Store store;
          for(; _first != _last; ++i_func, ++_first, ++i_coef, ++i_val )
          {
            T_COEF val = (*i_coef) * ( (*i_func)( *i_arg ) ) ;
            T_COEF grad = (*i_coef) * i_func->gradient( *i_arg );
            store.push_back( t_Pair( val, grad ) );
          }
          t_Store :: const_iterator i_val_begin = store.begin();
          t_Store :: const_iterator i_val_end = store.end();
          t_Store :: const_iterator i_val( i_val_begin );
          T_COEF result( 0 );
          for( types::t_unsigned N( store.size() ); N > ; N-- ) 
          {
            t_Store::const_iterator i_val( i_val_begin )
            types::t_unsigned n( store.size() )
            T_COEF intermediate( N != n ? i_val_begin->first: i_val_begin->second );
            for( --n, ++i_val; i_val ! = i_val_end; --n, ++i_val )
              result *= ( N != n ? i_val_begin->first: i_val_begin->second );
          }
          return result;
        }
      template< class T_COEF, class T_FUNC, class T_ARG, class T_RET>
        void operator( T_COEF _first, T_COEF _last, 
                       T_FUNC _func, T_ARG _arg, T_RET _ret ) 
        {
          // Stores evaluations and gradients.
          typedef std::pair< T_COEF, T_COEF > t_Pair;
          typedef std::vector< t_Pair > t_Store;
          t_Store store;
          for(; _first != _last; ++i_func, ++_first, ++i_coef, ++i_val )
          {
            T_COEF val = (*i_coef) * ( (*i_func)( *i_arg ) ) ;
            T_COEF grad = (*i_coef) * i_func->gradient( *i_arg );
            store.push_back( t_Pair( val, grad ) );
          }
          t_Store :: const_iterator i_val_begin = store.begin();
          t_Store :: const_iterator i_val_end = store.end();
          t_Store :: const_iterator i_val( i_val_begin );
          T_COEF result( 0 );
          for( types::t_unsigned N( store.size() ); N > ; N--, ++_ret) 
          {
            t_Store::const_iterator i_val( i_val_begin )
            types::t_unsigned n( store.size() )
            *_ret = ( N != n ? i_val_begin->first: i_val_begin->second );
            for( --n, ++i_val; i_val ! = i_val_end; --n, ++i_val )
              (*_ret) *= ( N != n ? i_val_begin->first: i_val_begin->second );
          }
        }
    };


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
       }



  } // end namespace details
  //! \endcond


}

#endif

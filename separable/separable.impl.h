//
//  Version: $Id$
//

namespace Separable
{
#define _FUNCBODY_ (code) \
        {\
          typedef typename t_Basis :: iterator function_iterator;\
          typedef typename t_Coefs :: const_iterator coef_iterator;\
          __ASSERT(     basis.size() == coefs.size()\
                    "Incoherent basis/coefficient sizes.\n" )\
          function_iterator i_func = basis.begin();\
          coef_iterator i_coef = coefs.begin();\
          t_Return result( (*i_coef) * ( (*i_func)( code _first ) ) ) ; \
          for( ++_first, ++i_func, ++i_coef; \
               _first != _last; ++i_func, ++_first, ++i_coef ) \
            groupop( result, scalarop( *i_coef, (*i_func)( code _first ) ) ); \
          return result; \
        }

  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITERATOR >
      typename Base<T_BASIS, T_LAW >::t_Return
        Base<T_BASIS, T_LAW> :: operator()( T_ITERATOR _first,
                                            T_ITERATOR _last ) const
        _FUNCBODY_( * )

  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITIN, class T_ITOUT>
      typename Base<T_BASIS, T_LAW> :: t_Return 
        Basis<T_BASIS, T_LAW> expand( T_ITIN _first, T_ITIN _last, T_ITOUT _out ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() == coefs.size()\
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( *_first ) ) ) ; 
          for( ++_first, ++i_func, ++i_coef; 
               _first != _last; ++i_func, ++_first, ++i_coef, ++_out ) 
            groupop( *_out, scalarop( *i_coef, (*i_func)( *_first ) ) ); 
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
      void BASE<T_BASIS, T_LAW> :: serialize( ARCHIVE & _ar,
                                              const unsigned int _version)
      {
        ar & basis;
        ar & coefs;
      }

  template< class T_BASIS, template<class> class T_LAW >
    template< class T_ITERATOR >
      typename Function<T_BASIS, T_LAW >::t_Return
        Function<T_BASIS, T_LAW> :: operator()( T_ITERATOR _first,
                                                T_ITERATOR _last ) const
        _FUNCBODY_()
        
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


  template< class T_ALLSQ, class T_BASIS >
    void assign_Allsq_matrix( const Function<T_BASIS> &_func,
                              const T_ALLSQ::t_Vectors &_input,
                              T_ALLSQ::t_Matrices &_A )
    {
      typedef T_ALLSQ t_Allsq;
      typedef typename t_Allsq :: t_Vectors t_Vectors;
      typedef typename t_Vectors :: t_Vector t_Vector;
      typedef typename t_Allsq :: t_Matrices t_Matrices;
      typedef const Function<T_BASIS> t_Function;

      // Variable expanded holder.
      typedef std::vector< t_Function::t_Basis::value_type::t_Basis::t_Return >  t_Vex;
      // Rank/Variable expanded holder.
      typedef std::vector< t_Vex >  t_Rex;
      // Observarion/Rank/Variable expanded holder.
      typedef std::vector< t_Rex > t_Obsex;
      // obsex will contain all factors multiplied but the collapsed one.
      // leading dimension is observations, followed by rank, and collapsed-dimension.
      t_Obsex obsex( input.size() ); 
      t_Obsex :: iterator i_ob = obsex.begin();
      t_Vectors :: const_iterator i_input = input.begin();
      t_Vectors :: const_iterator i_input_end = input.end();
      for(; i_input != i_input_end; ++i_input, ++i_ob )
      {
        i_ob->resize( _func.basis.size() );
        t_Rex :: iterator i_rank = i_ob->begin();
        t_Rex :: iterator i_rank_end = i_ob->end();
        t_Function::t_Basis::const_iterator i_funcrank = _func.basis.begin();
        for(; i_rank != i_rank_end; ++i_rank, ++i_func_rank )
        {
          t_Vex vex( i_funcrank->basis.size(), t_Vex::value_type(0) );
          i_funcrank->expand( i_input.begin(), i_input.end(), vex.begin() );

          t_Vex :: const_iterator i_vex = vex.begin();
          t_Vex :: const_iterator i_vex_end = vex.end();
          t_Vex :: iterator i_result = i_rank->begin();
          i_rank->resize( vex.size(), t_Vex::value_type(1) );
          // Multiplies all factors but one.
          for(; i_vex != i_vex_end; ++i_vex, ++i_funcfac )
          {
            t_Vex :: const_iterator i_vex2 = vex.begin();
            for(; i_vex2 != i_vex_end; ++i_vex2 )
              if( i_vex2 != i_vex ) *i_result *= i_vex2;
          }
        }
      }
        
      // Organization of the _A matrices.
      // Each dimension dependent sub-matrix is one value in the vector _A.
      // each sub-matrix is a vector organized as follows
      //   obs 0     obs 1
      // yyyyyyyy  yyyyyyyyyy
      // Each yyyyyyy runs as 
      //   r=0    r=1  rank of separable function
      //  xxxxx  xxxxx 
      //  Each xxxxx as 
      //   m=0     m=1   basis for factor n of separable function.
      //  aaaaaa aaaaaa
      t_Allsq :: t_Matrices :: iterator i_mat = _A.begin();
      t_Allsq :: t_Matrices :: iterator i_mat_end = _A.end();
      // Loop over collapsed dimensions.
      for(types::t_int dim=0; i_mat != i_mat_end; ++i_mat, ++dim )
      {
        // Loop over observations
        i_input = i_input.begin();
        i_ob = obsex.begin();
        i_mat->reserve( input.size() *  obsex.size()
                                     *  obsex[0].size()
                                     *  obsex[0][0].size() );
        for(; i_input != i_input_end; ++i_input, ++obsex )
        {


          // Loop over ranks. 
          t_Function::t_Basis::const_iterator i_rank = _func.basis.begin();
          t_Function::t_Basis::const_iterator i_rank_end = _func.basis.end();
          t_Rex :: iterator i_rex = i_ob->begin();
          for(; i_rank != i_rank_end; ++i_rank, ++i_rex )
          {
            // Finds function for this dimension
            const t_Function::t_Basis::value_type &dim_func = (*i_rank)[dim];
            // Finds factor for this dimension
            t_Vex :: value_type factor( (*i_rex)[dim] );

            // Loops over 1d basis of this collapsed dimension.
            t_Function::t_Basis::value_type::t_Basis::const_iterator i_1d = dim_func->begin();
            t_Vector :: const_iterator i_var = i_input->begin();
            t_Vector :: const_iterator i_var_end = i_input->begin();
            for(; i_var != i_var_end; ++i_1d, ++i_var )
              i_mat->push_back( (*i_1d)( *i_var ) * factor );  

          } // end loop over ranks.

        }  // end loop over observations

      } // end loop over collapsed dimensions
    }

  template< class T_ALLSQ, class T_BASIS >
    assign_from_allsq<T_BASIS,T_ALLSQ>( Function<T_BASIS> &_func,
                                        const T_ALLSQ::t_Vectors &_coefs )
    {
      // Type of the fitting method.
      typedef T_ALLSQ t_Allsq;
      // Type of the input vectors.
      typedef typename t_Allsq :: t_Vectors t_Vectors;
      // Type of a single input vector.
      typedef typename t_Vectors :: t_Vector t_Vector;

      t_Vectors :: const_iterator i_vec = _coefs.begin();
      t_Vectors :: const_iterator i_vec_end = _coefs.end();
      // Loop over collapsed dimensions.
      for(types::t_int dim=0; i_vec != i_vec_end; ++i_vec, ++dim )
      {
        // Loop over ranks. 
        t_Function::t_Basis::const_iterator i_rank = _func.basis.begin();
        t_Function::t_Basis::const_iterator i_rank_end = _func.basis.end();
        for(; i_rank != i_rank_end; ++i_rank )
        {
          // Finds function for this dimension
          const t_Function::t_Basis::value_type &dim_func = (*i_rank)[dim];
         
          // Loops over 1d basis of this collapsed dimension.
          t_Vector :: const_iterator i_var = i_vec->begin();
          t_Vector :: const_iterator i_var_end = i_vec->begin();
          for(types::t_int i=0; i_var != i_var_end; ++i, ++i_var )
            dim_func[i] = *i_var;

        } // end loop over ranks.

      } // end loop over collapsed dimensions
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

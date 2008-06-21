//
//  Version: $Id$
//

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
          __ASSERT(     basis.size() == coefs.size(),
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( *_first ) ) ) ; 
          for( ++_first, ++i_func, ++i_coef; 
               _first != _last; ++i_func, ++_first, ++i_coef ) 
            groupop( result, scalarop( *i_coef, (*i_func)( *_first ) ) ); 
          return result; 
        }

  template< class T_B, template<class> class T_G, template<class> class T_S >
    template< class T_ITIN, class T_ITOUT>
      typename Base<T_B,T_G,T_S> :: t_Return  
        Base<T_B,T_G,T_S> :: expand( T_ITIN _first, T_ITIN _last, T_ITOUT _out ) const
        {
          typedef typename t_Basis :: iterator function_iterator;
          typedef typename t_Coefs :: const_iterator coef_iterator;
          __ASSERT(     basis.size() == coefs.size(),
                    "Incoherent basis/coefficient sizes.\n" )
          function_iterator i_func = basis.begin();
          coef_iterator i_coef = coefs.begin();
          t_Return result( (*i_coef) * ( (*i_func)( *_first ) ) ) ; 
          for( ++_first, ++i_func, ++i_coef; 
               _first != _last; ++i_func, ++_first, ++i_coef, ++_out ) 
            groupop( *_out, scalarop( *i_coef, (*i_func)( *_first ) ) ); 
          return result; 
        }

  template< class T_B, template<class> class T_G, template<class> class T_S >
    typename Base<T_B,T_G,T_S>::t_Return
      Base<T_B,T_G,T_S> :: operator()( t_Arg _arg ) const
      {
        typedef typename t_Basis :: const_iterator function_iterator;
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

  template< class T_B, template<class> class T_G, template<class> class T_S >
    template< class T_ITERATOR >
      void Base<T_B,T_G,T_S> :: set( T_ITERATOR _first, T_ITERATOR _last )
      {
        __ASSERT( (_last-_first) == coefs.size(), 
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
        __ASSERT(     basis.size() == coefs.size(),
                  "Incoherent basis/coefficient sizes.\n" )
        function_iterator i_func = basis.begin();
        function_iterator i_func_end = basis.end();
        coef_iterator i_coef = coefs.begin();
        t_Return result( (*i_coef) * ( (*i_func)( _first, _last ) ) ) ; 
        for( ++i_func, ++i_coef; i_func != i_func_end; ++i_func, ++i_coef ) 
          groupop( result, scalarop( *i_coef, (*i_func)( _first, _last ) ) ); 
        return result; 
      }
        
  template<class T_FUNCTION> template< class T_VECTORS >
  void Collapse<T_FUNCTION> :: init( const T_VECTORS &_x ) 
  {
    // finds maximum size of separable function.
    {
      typename t_Function::t_Basis::const_iterator i_first = function.basis.begin();
      typename t_Function::t_Basis::const_iterator i_last = function.basis.end();
      for(D=0; i_first != i_last; ++i_first )
        if( D < i_first->basis.size() ) D = i_first->basis.size();
    }

    // computes sizes
    sizes.resize( D );
    for( size_t d = 0; d < D; ++d )
    {
      typedef typename t_Function::t_Basis::const_iterator t_const_iterator;
      t_const_iterator i_rank = function.basis.begin();
      t_const_iterator i_rank_end = function.basis.end();
      for(; i_rank != i_rank_end; ++i_rank )
      {
        if( i_rank->basis.size() < d ) continue;
        sizes[d].push_back( i_rank->basis[d].basis.size() );
      }
    }
                      
    // computes expanded matrix.
    nb_obs = _x.size();
    nb_ranks = function.basis.size();
    expanded.reserve( nb_obs * D );
    typename T_VECTORS::const_iterator i_x = _x.begin(); 
    typename T_VECTORS::const_iterator i_x_end = _x.end(); 
    for(; i_x != i_x_end; ++i_x ) // loop over observables.
    {
      typename T_VECTORS::value_type::const_iterator i_var = i_x->begin();
      for( size_t d = 0; d < D; ++d, ++i_var )
      {
        __ASSERT( i_var != i_x->end(), "Inconsistent sizes.\n" )

        expanded.resize( expanded.size() + 1 );

        typename t_Expanded :: value_type &row = expanded.back();
        typedef typename t_Function::t_Basis::const_iterator t_const_iterator;
        t_const_iterator i_rank = function.basis.begin();
        t_const_iterator i_rank_end = function.basis.end();
        for(; i_rank != i_rank_end; ++i_rank ) // loop over ranks
        {
          if( i_rank->basis.size() < d ) continue;
          typedef typename t_Function::t_Basis::value_type
                                     ::t_Basis::const_iterator t_citerator;
          t_citerator i_func = i_rank->basis.begin();
          t_citerator i_func_end = i_rank->basis.end();
          for(; i_func != i_func_end; ++i_func )
            row.push_back( (*i_func)(*i_var) );
        } // end of loop over ranks.
      } // end of loop over dimensions.
    } // end of loop over observables.
  }

  template<class T_FUNCTION>  template< class T_MATRIX, class T_VECTORS >
    void Collapse<T_FUNCTION>::operator()( T_MATRIX &_A, types::t_unsigned _dim,
                                           const T_VECTORS &_coefs )
    {
      if( ( not is_initialized ) or ( ( not do_update ) and _dim == 0 ) )
        initialize_factors( _coefs );
      else if ( do_update ) update_factors( _dim, _coefs );
      is_initialized = true;   
 
      _A.resize( nb_obs );
      typename T_MATRIX :: iterator i_A = _A.begin();
      typename T_MATRIX :: iterator i_A_end = _A.end();
      typename t_Expanded :: const_iterator i_ex = expanded.begin() + _dim;
      typename t_Factors :: const_iterator i_facs = factors.begin();
      for(; i_A != i_A_end; i_ex += _dim ) // loops over observables
      {
        t_Sizes :: value_type :: const_iterator i_size = sizes[_dim].begin();
        t_Sizes :: value_type :: const_iterator i_size_end = sizes[_dim].end();
        __ASSERT( factors.size() != sizes[_dim].size(),
                  "Inconsistent sizes.\n" )
        typename t_Expanded :: value_type :: const_iterator i_var = i_ex->begin();
        for(; i_size != i_size_end; ++i_size, ++i_facs ) // loops over ranks.
        {
          // Computes factor for this rank.
          typename t_Factors :: value_type :: const_iterator i_fac = i_facs->begin();
          typename t_Factors :: value_type :: const_iterator i_fac_end = i_facs->end();
          typename t_Factors ::value_type :: value_type factor(1);
          for(types::t_unsigned i = 0; i_fac != i_fac_end; ++i_fac, ++i )
            if( i != _dim ) factor *= (*i_fac);
 
          // loops over basis set of factor function of current dimension and rank.
          for( types::t_unsigned i = *i_size; i > 0; --i, ++i_var, ++i_A )
          {
            __ASSERT( i_var != i_ex->end(), "Inconsistent sizes.\n" )
            __ASSERT( i_A != _A.end(), "Inconsistent sizes.\n" )
            *i_A = (*i_var) * factor;
          }
        }
      }
    }

  template< class T_FUNCTION> template< class T_VECTORS >
    void Collapse<T_FUNCTION>::initialize_factors( const T_VECTORS &_coefs )
    {
      __ASSERT( nb_ranks != function.basis.size(),
                "Inconsistent rank size.\n" )
      factors.resize( nb_ranks * nb_obs );
      typename t_Function :: t_Basis :: const_iterator i_rank=function.basis.begin();
      typename t_Factors :: iterator i_facs = factors.begin();
      typename t_Factors :: iterator i_facs_end = factors.end();
      for(types :: t_unsigned o(0); o < nb_obs; --o ) // over observables
      {
        // loop over ranks
        for(types::t_unsigned r(0); i_facs != i_facs_end; ++i_rank, ++i_facs, ++r ) 
        {
          i_facs->resize( i_rank->basis.size() );
          // loop over dimensions
          typename t_Factors :: value_type :: iterator i_fac = i_facs->begin();
          typename t_Factors :: value_type :: iterator i_fac_end = i_facs->end();
          for(types::t_unsigned d = 0; i_fac != i_fac_end; ++i_fac, ++d )
          {
            *i_fac = t_Type(0);
            // performs cdot calculation: 
            //  _ first compute the accumulated number of elements over ranks
            //    and basis functions.
            __ASSERT( sizes.size() <= d,
                      "Inconsistent size of array Collapse::sizes.\n" )
            __ASSERT( sizes[d].size() <= r,
                      "Inconsistent size of array Collapse::sizes.\n" )
            types::t_unsigned acc = std::accumulate( sizes[d].begin(),
                                                     sizes[d].begin() + r,
                                                     types::t_unsigned(0) );
            //  _ second retrieve iterator to expanded function for current
            //    observable, rank, and dimension.
            __ASSERT( expanded.size() <= d + o * D,
                      "Inconsistent size of input vectors.\n" )
            __ASSERT( expanded[d+o*D].size() < acc + sizes[d][r],
                      "Inconsistent size of input vectors.\n" )
            typedef typename t_Expanded :: value_type :: const_iterator t_citerator;
            t_citerator i_func =   expanded[ d + o * D ].begin()  + acc;
            t_citerator i_func_end =  i_func + sizes[d][r];
            // _ second retrieve the coefficients.
            typedef typename T_VECTORS :: value_type :: const_iterator t_coefit; 
            __ASSERT( _coefs.size() <= d, "Inconsistent size of input vectors." )
            __ASSERT( _coefs[d].size() < acc+sizes[d][r],
                      "Inconsistent size of input vectors." )
            t_coefit i_coef = _coefs[d].begin() + acc;
            // _ finally perform dot product and store in factors matrix.
            for(; i_func != i_func_end; ++i_func, i_coef )
              *i_fac += ( *i_coef  ) * ( *i_func );
          } // end of loop over dimensions
        } // end of loop over ranks
      } // end of loop over observables
    }
 
  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION>::update_factors( types::t_unsigned _dim,
                                               const T_VECTORS &_coefs )
    {
      _dim == 0 ? _dim = D: --_dim;
      __ASSERT( nb_ranks != function.basis.size(),
                "Inconsistent rank size.\n" )
      factors.resize( nb_ranks * nb_obs );
      typename t_Function :: t_Basis :: const_iterator i_rank=function.basis.begin();
      typename t_Factors :: iterator i_facs = factors.begin();
      typename t_Factors :: iterator i_facs_end = factors.end();
      for(types :: t_unsigned o(0); o < nb_obs; --o ) // over observables
      {
        // loop over ranks
        for(types::t_unsigned r(0); i_facs != i_facs_end; ++i_rank, ++i_facs, ++r ) 
        {
          __ASSERT( i_facs->size() <= _dim, "Incoherent number of factors.\n" )
          t_Type &factor = (*i_facs)[_dim];
          factor = t_Type(0);
          // performs cdot calculation: 
          //  _ first compute the accumulated number of elements over ranks and
          //    basis functions.
          __ASSERT( sizes.size() <= _dim,
                    "Inconsistent size of array Collapse::sizes.\n" )
          __ASSERT( sizes[_dim].size() <= r,
                    "Inconsistent size of array Collapse::sizes.\n" )
          types::t_unsigned acc = std::accumulate( sizes[_dim].begin(),
                                                   sizes[_dim].begin() + r,
                                                   types::t_unsigned(0) );
          //  _ second retrieve iterator to expanded function for current
          //    observable, rank, and dimension.
          __ASSERT( expanded.size() <= _dim + o * D,
                    "Inconsistent size of input vectors.\n" )
          __ASSERT( expanded[_dim+o*D].size() < acc + sizes[_dim][r],
                    "Inconsistent size of input vectors.\n" )
          typedef typename t_Expanded :: value_type :: const_iterator t_citerator;
          t_citerator i_func =   expanded[ _dim + o * D ].begin()  + acc;
          t_citerator i_func_end =  i_func + sizes[_dim][r];
          // _ second retrieve the coefficients.
          typedef typename T_VECTORS :: value_type :: const_iterator t_coefit; 
          __ASSERT( _coefs.size() <= _dim, "Inconsistent size of input vectors." )
          __ASSERT( _coefs[_dim].size() < acc+sizes[_dim][r],
                    "Inconsistent size of input vectors." )
          t_coefit i_coef = _coefs[_dim].begin() + acc;
          // _ finally perform dot product and store in factors matrix.
          for(; i_func != i_func_end; ++i_func, i_coef )
            factor += ( *i_coef  ) * ( *i_func );
        } // end of loop over ranks
      } // end of loop over observables
    }

  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION> :: reassign( const T_VECTORS& _solution ) const
    {
      typedef T_VECTORS t_Vectors;
      typedef typename t_Vectors :: const_iterator const_solit;
      typedef typename t_Vectors :: value_type :: const_iterator const_soldimit;
      const_solit i_dim = _solution.begin();
      const_solit i_dim_end = _solution.end();
      for(size_t d=0; i_dim != i_dim_end; ++i_dim, ++d ) // loop over dimensions.
      {
        const_soldimit i_coef = i_dim->begin();
        const_soldimit i_coef_ = i_coef;
        const_soldimit i_coef_end = i_dim->begin();
        for(size_t r=0; i_coef != i_coef_end; ++r )
        {
          typedef typename t_Function :: t_Basis :: value_type
                                      :: t_Basis :: value_type t_factor_func;
          t_factor_func &facfunc = function.basis[r].basis[d];
          i_coef_ += facfunc.coefs.size();
          __ASSERT( i_coef_end - i_coef_ < 0,
                    "Inconsistent number of coefficients." )
          std::copy( i_coef, i_coef_, facfunc.coefs.begin() );
          i_coef = i_coef_;
        }
      }
    }

  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION> :: create_coefs( T_VECTORS& _coefs ) const
    {
      typedef T_VECTORS t_Vectors;
      typedef typename t_Vectors :: iterator t_solit;
      _coefs.resize( sizes.size() );
      t_solit i_dim = _coefs.begin();
      t_solit i_dim_end = _coefs.end();
      t_Sizes :: const_iterator i_size = sizes.begin();
      for(; i_dim != i_dim_end; ++i_dim, ++i_size ) // loop over dimensions.
      {
        types::t_int ri = std::accumulate( i_size->begin(), i_size->end(), 0 );
        i_dim->resize( ri );
      }
    }

  template< class T_FUNCTION >
    void Collapse<T_FUNCTION> :: reset() 
    {
      is_initialized = false;
      expanded.clear();
      factors.clear();
      sizes.clear();
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

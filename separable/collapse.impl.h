//
//  Version: $Id$
//


namespace Separable
{
  template<class T_FUNCTION> template< class T_VECTORS >
  void Collapse<T_FUNCTION> :: init( const T_VECTORS &_x ) 
  {
    // finds maximum size of separable function.
    typename t_Function::t_Basis::const_iterator i_first = function.basis.begin();
    typename t_Function::t_Basis::const_iterator i_last = function.basis.end();
    for(D=0; i_first != i_last; ++i_first )
      if( D < i_first->basis.size() ) D = i_first->basis.size();

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
        __ASSERT( i_var == i_x->end(), "Inconsistent sizes.\n" )

        expanded.resize( expanded.size() + 1 );

        typename t_Expanded :: value_type &row = expanded.back();
        typedef typename t_Function::t_Basis::const_iterator t_const_iterator;
        t_const_iterator i_rank = function.basis.begin();
        t_const_iterator i_rank_end = function.basis.end();
        for(; i_rank != i_rank_end; ++i_rank ) // loop over ranks
        {
          if( i_rank->basis.size() < d ) continue;
          typedef typename t_Function::t_Basis::value_type
                                     ::t_Basis::value_type
                                     ::t_Basis::const_iterator t_citerator;
          t_citerator i_func = i_rank->basis[d].basis.begin();
          t_citerator i_func_end = i_rank->basis[d].basis.end();
          for(; i_func != i_func_end; ++i_func )
            row.push_back( (*i_func)(*i_var) );
        } // end of loop over ranks.
      } // end of loop over dimensions.
    } // end of loop over observables.
  }

#ifdef _LADADEBUG
  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION> :: print_funcs( const T_VECTORS &_x ) const
    {
      typename T_VECTORS::const_iterator i_x = _x.begin(); 
      typename T_VECTORS::const_iterator i_x_end = _x.end(); 
      for(; i_x != i_x_end; ++i_x ) // loop over observables.
      {
        typename T_VECTORS::value_type::const_iterator i_var = i_x->begin();
        std::cout << "Funcs: " << std::endl;
        typedef typename t_Function::t_Basis::const_iterator t_rank; 
        t_rank i_rank = function.basis.begin();
        t_rank i_rank_end = function.basis.end();
        for(; i_rank != i_rank_end; ++i_rank )
        {
          typedef typename t_rank::value_type::t_Basis::const_iterator t_dim; 
          t_dim i_dim = i_rank->basis.begin();
          t_dim i_dim_end = i_rank->basis.end();
          typename T_VECTORS::value_type::const_iterator i_var = i_x->begin();
          for(; i_dim != i_dim_end; ++i_dim, ++i_var )
          {
            __ASSERT( i_var == i_x->end(), "Inconsistent sizes.\n" )
            using namespace boost::lambda;
            std::for_each
            ( 
              i_dim->basis.begin(), i_dim->basis.end(),
              std::cout << bind<types::t_real>( _1, constant(*i_var) ) << " "
            );
          } // dim
          std::cout << "\n";
        } // rank
        std::cout << "\n\n";
      }
    }
#endif

  template<class T_FUNCTION>  template< class T_MATRIX, class T_VECTORS >
    void Collapse<T_FUNCTION>::operator()( T_MATRIX &_A, types::t_unsigned _dim,
                                           const T_VECTORS &_coefs )
    {
#     ifdef _LADADEBUG
        try {
#     endif
      __ASSERT( sizes.size() <= _dim, "Inconsistent sizes.\n" )
      __ASSERT( _coefs.size() != D,
                "Inconsistent number of dimensions/coefficients.\n" )
      __ASSERT
      ( 
        _coefs[_dim].size() != std::accumulate( sizes[_dim].begin(),
                                                sizes[_dim].end(), 0 ),
        "Inconsistent number of ranks/dimensions/coefficients.\n" 
      )
      __ASSERT( factors.size() == _coefs[_dim].size(), "Inconsistent sizes.\n" )
      __ASSERT( expanded.size() <= _dim, "Inconsistent sizes.\n" )

      if( ( not is_initialized ) or ( ( not do_update ) and _dim == 0 ) )
        initialize_factors( _coefs );
      else if ( do_update ) update_factors( _dim, _coefs );
      is_initialized = true;   
 
      _A.resize( nb_obs * _coefs[_dim].size() );
      typename T_MATRIX :: iterator i_A = _A.begin();
      typename T_MATRIX :: iterator i_A_end = _A.end();
      typename t_Expanded :: const_iterator i_ex = expanded.begin() + _dim;
      typename t_Factors :: const_iterator i_facs = factors.begin();
      std::cout << "C0 A:\n";
      for(; i_A != i_A_end; i_ex += D ) // loops over observables
      {
        __ASSERT( i_ex->size() != _coefs[_dim].size(), "Inconsistent sizes.\n" )
        t_Sizes :: value_type :: const_iterator i_size = sizes[_dim].begin();
        t_Sizes :: value_type :: const_iterator i_size_end = sizes[_dim].end();
        typename t_Expanded :: value_type :: const_iterator i_var = i_ex->begin();
        for(; i_size != i_size_end; ++i_size, ++i_facs ) // loops over ranks.
        {
          // Computes factor for this rank.
          typename t_Factors :: value_type :: const_iterator i_fac = i_facs->begin();
          typename t_Factors :: value_type :: const_iterator i_fac_end = i_facs->end();
          typename t_Factors ::value_type :: value_type factor(1);
          for(types::t_unsigned d = 0; i_fac != i_fac_end; ++i_fac, ++d )
            if( d != _dim ) factor *= (*i_fac);
 
          // loops over basis set of factor function of current dimension and rank.
          for( types::t_unsigned i = *i_size; i > 0; --i, ++i_var, ++i_A )
            { *i_A = (*i_var) * factor; std::cout << *i_A << " "; }
        } // end of loop over ranks.
      } // end of loop over observables.
      std::cout << "\nC A: ";
      std::for_each( _A.begin(), _A.end(), std::cout << boost::lambda::_1 << " " );
      std::cout << "\n";
#     ifdef _LADADEBUG
        } __CATCHCODE(, "Error while creating collapsed A-matrix.\n" )
#     endif
    }

  template< class T_FUNCTION> template< class T_VECTORS >
    void Collapse<T_FUNCTION>::initialize_factors( const T_VECTORS &_coefs )
    {
#     ifdef _LADADEBUG
        try {
#     endif
      __ASSERT( nb_ranks != function.basis.size(),
                "Inconsistent rank size.\n" )
      factors.resize( nb_ranks * nb_obs );
      typename t_Factors :: iterator i_facs = factors.begin();
      typename t_Factors :: iterator i_facs_end = factors.end();
      for(types :: t_unsigned o(0); o < nb_obs; ++o ) // over observables
      {
        // loop over ranks
        typedef typename t_Function :: t_Basis :: const_iterator t_rankit;
        t_rankit i_rank=function.basis.begin();
        t_rankit i_rank_end=function.basis.end();
        for(types::t_unsigned r(0); i_rank != i_rank_end; ++i_rank, ++i_facs, ++r ) 
        {
          i_facs->resize( i_rank->basis.size(), t_Type(0) );
          // loop over dimensions
          typename t_Factors :: value_type :: iterator i_fac = i_facs->begin();
          typename t_Factors :: value_type :: iterator i_fac_end = i_facs->end();
          for(types::t_unsigned d = 0; i_fac != i_fac_end; ++i_fac, ++d )
          {
            // performs cdot calculation: 
            //  _ first compute the accumulated number of elements over ranks
            //    and basis functions.
            __ASSERT( sizes.size() <= d,
                      "Inconsistent size of array Collapse::sizes.\n" )
            __ASSERT( sizes[d].size() <= r,
                      "Inconsistent size of array Collapse::sizes[" << d << "].\n" )
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
            for(; i_func != i_func_end; ++i_func, ++i_coef )
              *i_fac += ( *i_coef  ) * ( *i_func );
          } // end of loop over dimensions
        } // end of loop over ranks
      } // end of loop over observables
#     ifdef _LADADEBUG
        } __CATCHCODE(, "Error while creating factors.\n" )
#     endif
    }
 
  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION>::update_factors( types::t_unsigned _dim,
                                               const T_VECTORS &_coefs )
    {
#     ifdef _LADADEBUG
        try {
#     endif
      _dim == 0 ? _dim = D-1: --_dim;
      __ASSERT( nb_ranks != function.basis.size(),
                "Inconsistent rank size.\n" )
      __ASSERT( factors.size() != nb_ranks * nb_obs,
                "Inconsistent factors size.\n" )
      typename t_Factors :: iterator i_facs = factors.begin();
      typename t_Factors :: iterator i_facs_end = factors.end();
      for(types :: t_unsigned o(0); o < nb_obs; ++o ) // over observables
      {
        // loop over ranks
        typedef typename t_Function :: t_Basis :: const_iterator t_rankit;
        t_rankit i_rank=function.basis.begin();
        t_rankit i_rank_end=function.basis.end();
        for(types::t_unsigned r(0); i_rank != i_rank_end; ++i_rank, ++i_facs, ++r ) 
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
          for(; i_func != i_func_end; ++i_func, ++i_coef )
            factor += ( *i_coef  ) * ( *i_func );
        } // end of loop over ranks
      } // end of loop over observables
#     ifdef _LADADEBUG
        } __CATCHCODE(, "Error while updating factors.\n" )
#     endif
    }

  template< class T_FUNCTION > template< class T_VECTORS >
    void Collapse<T_FUNCTION> :: reassign( const T_VECTORS& _solution ) const
    {
#     ifdef _LADADEBUG
        try {
#     endif
      typedef T_VECTORS t_Vectors;
      typedef typename t_Vectors :: const_iterator const_solit;
      typedef typename t_Vectors :: value_type :: const_iterator const_soldimit;
      const_solit i_dim = _solution.begin();
      const_solit i_dim_end = _solution.end();
      for(size_t d=0; i_dim != i_dim_end; ++i_dim, ++d ) // loop over dimensions.
      {
        const_soldimit i_coef = i_dim->begin();
        const_soldimit i_coef_ = i_coef;
        const_soldimit i_coef_end = i_dim->end();
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
#     ifdef _LADADEBUG
        } __CATCHCODE(, "Error while assigning solution coefs to function.\n" )
#     endif
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

}

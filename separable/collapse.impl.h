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
    // finds nb of target values
    nb_targets = _x.size();
    // finds nb of ranks
    nb_ranks = function.basis.size();

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

    // Computes expanded
    expanded.resize( D );
    typename t_Expanded :: iterator i_dexp = expanded.begin();
    t_Sizes :: const_iterator i_dsize = sizes.begin();
    for( types::t_unsigned d(0); d < D; ++d, ++i_dexp, ++i_dsize, ++i_dfacs )
    {
      // prepares loop over ranks
      i_dexp->resize( rsize.end() );
      typename t_Expanded :: value_type :: iterator i_rexp = i_dexp->begin();
      typename t_Sizes :: value_type :: const_iterator i_rsize = i_dsize->begin();
      typename t_Sizes :: value_type :: const_iterator i_rsize_end = i_dsize->end();
      for( types::t_unsigned r(0); i_rsize != i_rsize_end ; ++r, ++i_rexp, ++i_rsize ) 
      {
        __ASSERT( _function.basis[r].basis[d].size() != *i_rsize,
                  "Inconsistent size.\n" )

        // prepares loop over functions.
        i_rex->resize( (*i_rsize) * nb_obs );
        typename t_Expanded :: value_type
                            :: value_type iterator i_iexp = i_rexp->begin();
        typename t_Function :: t_Basis :: value_type
                            :: t_Basis :: value_type 
                            :: t_Basis :: const_iterator i_ifunc = _function.basis[r]
                                                                            .basis[d].begin();
        namespace bl = boost::lambda;
        for( types::t_unsigned i(0); i < *i_rsize; ++i, ++i_ifunc )
          opt::concurrent_loop // loop over target values
          (
            _x.begin(), _x.end(), boost::ref( i_iexp ),
            // computes expanded[d,r, (i,o)] = g_{d,i}^{(r)}(x_d^{(o)}).
            bl::_2 = bl::bind( bl::var( *i_ifunc ), bl::_1 )
          ); // end of loop over target values


        __ASSERT( i_iexp  != i_rexp->end(), "Inconsistent size.\n" )
        __ASSERT( i_ifunc != i_function.basis[r].basis[d].end(), "Inconsistent size.\n" )
      } // end of loop over ranks
      __ASSERT( i_rexp  != i_dexp->end(), "Inconsistent size.\n" )
      __ASSERT( i_rsize != i_dsize->end(), "Inconsistent size.\n" )


    } // end of loop over dimensions.
    __ASSERT( i_dexp  != expanded.end(), "Inconsistent size.\n" )
    __ASSERT( i_dsize != sizes.end(), "Inconsistent size.\n" )



    // computes factors.
  }

  template<class T_FUNCTION>  template< class T_MATRIX, class T_VECTORS >
    void Collapse<T_FUNCTION>::matrixA( T_MATRIX &_A, types::t_unsigned _dim,
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
      for( types::t_unsigned r(0); r < nb_ranks; ++r )
      {
        // loop over target values.
        for(types::t_unsigned o(0); o < nb_targets; ++o, ++i_facs )
        {
          i_facs->resize( D );
          typename t_Factors :: value_type :: iterator i_fac = i_facs->begin();
          for( types::t_unsigned d(0); d < D; ++d, ++i_fac )
          {
            __ASSERT( function.size() <= r, "Inconsistent rank size.\n" )
            if( function.basis[r].size() <= d ) { *i_fac = t_Type(1); continue; }

            types::t_unsigned acc = std::accumulate( sizes[d].begin(),
                                                     sizes[d].begin() + r, 0 );
            __ASSERT( _coefs.size()  <= d, "Inconsistent dimension size.\n" )
            __ASSERT( _coefs[d].size() <= acc + sizes[d][r], "Inconsistent size.\n" )
            __ASSERT( expanded[d][r].size() <= ( acc + sizes[d][r] ) * nb_targets + o, 
                       "Inconsistent size.\n" )

            typedef typename t_Expanded :: value_type :: value_type :: const_iterator t_itexp;
            typedef typename T_VECTORS :: value_type :: const_iterator t_itcoefs;
     
            t_coefs i_coef = _coefs[d].begin() + acc;
            t_itexp i_iexp = expanded[d][r].begin() + acc * nb_targets + o;

            *i_fac = t_Type(0);
            for( types::t_unsigned i( sizes[d][r] ); i > o; ++i, i_iexp += o, ++i_coef )
              *i_fac += (*i_coef) * (*i_iexp ); 
            
          } // end of loop over dimensions
        } // end of loop over target values.
        __ASSERT( i_fac != factors.end(), "Inconsistent size.\n" )
      } // end of loop over ranks
      __ASSERT( i_facs != factors.end(), "Inconsistent size.\n" )

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

      __ASSERT( nb_ranks != function.basis.size(),
                "Inconsistent rank size.\n" )
      __ASSERT( factors.size() != nb_ranks * nb_obs,
                "Inconsistent rank size.\n" )

      _dim == 0 ? _dim = D-1: --_dim;
      typename t_Factors :: iterator i_facs = factors.begin();
      for( types::t_unsigned r(0); r < nb_ranks; ++r )
      {
        // loop over target values.
        for(types::t_unsigned o(0); o < nb_targets; ++o, ++i_facs )
        {
          __ASSERT( factors.size() != D, "Inconsistent size.\n" )
          if( function.basis[r].size() <= _dim ) { (*i_facs)[_dim] = t_Type(1); continue; }

          types::t_unsigned acc = std::accumulate( sizes[_dim].begin(),
                                                   sizes[_dim].begin() + r, 0 );
          __ASSERT( _coefs.size()  <= _dim, "Inconsistent dimension size.\n" )
          __ASSERT( _coefs[d].size() <= acc + sizes[_dim][r], "Inconsistent size.\n" )
          __ASSERT( expanded[d][r].size() <= ( acc + sizes[_dim][r] ) * nb_targets + o, 
                     "Inconsistent size.\n" )

          typedef typename t_Expanded :: value_type :: value_type :: const_iterator t_itexp;
          typedef typename T_VECTORS :: value_type :: const_iterator t_itcoefs;
     
          t_coefs i_coef = _coefs[_dim].begin() + acc;
          t_itexp i_iexp = expanded[_dim][r].begin() + acc * nb_targets + o;

          (*i_facs)[_dim] = t_Type(0);
          for( types::t_unsigned i( sizes[_dim][r] ); i > o; ++i, i_iexp += o, ++i_coef )
            (*i_facs)[_dim] += (*i_coef) * (*i_iexp ); 
          
        } // end of loop over target values.
      } // end of loop over ranks
      __ASSERT( i_facs != factors.end(), "Inconsistent size.\n" )

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

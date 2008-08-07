
namespace CE
{

  template< class T_SOLVER >
    opt::ErrorTuple Fit :: operator( t_Vector &_x, T_SOLVER &_solver ) const
    {
      __DEBUGTRYBEGIN
      namespace bl = boost::lambda;
      namespace bblas = boost::numeric::ublas;
      __ASSERT( nb_cls != clusters.size(),
                "Inconsistent number of clusters.\n" )
      __ASSERT( targets.size() != weights.size(),
                "Inconsistent number of targets and weights.\n" )
      __ASSERT( pis.size() != targets.size(),
                "Inconsistent number of targets and pis.\n" )
 
      t_Matrix A( nb_cls, nb_cls );
      t_Vector b( nb_cls );
      create_A_n_b( A, b );
      solver( A, _x, b );
 
      // computes square errors.
      return check_training( _x, verbose );
      __DEBUGTRYEND(, "Error in Fit::operator().\n" )
    }

  template< class T_POLICY >
  void Fit<T_POLICY> :: add_to_A_n_b( t_Vector &_A, t_Vector &_b,
                                      const t_StructPis &_pis,
                                      const types::t_real _weight,
                                      const types::t_real _energy )
  {
    __ASSERT( i_pis->size() != nb_cls, "Inconsistent number of targets pis.\n" )
    // loop over alpha.
    t_ESums :: iterator i_b = _b.begin();
    t_PSums :: array_type :: iterator i_A = _A.data().begin();
    t_StructPis :: const_iterator i_pi_begin = _pis.begin();
    t_StructPis :: const_iterator i_pi_end = i_pis.end();
    t_StructPis :: const_iterator i_alphapi( i_pi_begin );
    for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_esum)
    {
      __ASSERT( i_esum == _b.end(), "Iterator out of range.\n" )
      *i_b += _weight * energy * (*i_alphapi);

      // loop over betas.
      t_StructPis :: const_iterator i_betapi( i_pi_begin );
      for(; i_betapi != i_pi_end; ++i_betapi, ++i_A )
      {
        __ASSERT( i_A == _A.data().end(), "Iterator out of range.\n" )
        *i_A += _weight * (*i_alphapi) * (*i_betapi);
      }
    } // end of loop over alphas.
  } 

  template<class T_POLICY>
  void Fit<T_POLICY> :: create_A_n_b( t_Vector &_A, t_Vector &_b )
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
    __ASSERT( pis.size() != structures.size(),
              "Inconsistent number of structures and pis.\n" )
    __ASSERT( structures.size() != weights.size(),
              "Inconsistent number of structures and weights.\n" )

    // Resizes the clusters and fills with zeros.
    _b.resize( nb_cls );
    _A.resize( nb_cls, nb_cls );
    std::fill( _b.begin(), _b.end(), 0e0 );
    std::fill( _A.data().begin(), _A.data().end(), 0e0 );

    // loop over targets.
    t_Structures :: const_iterator i_target = structures.begin();
    t_Structures :: const_iterator i_target_end = structures.end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    for( types::t_unsigned i(0); 
         i_target != i_target_end; 
         ++i_target, ++i_pis, ++i_w, ++i )
    {
      if( t_Policy :: found( i ) ) continue; 
      add_to_A_n_b( _A, _b, *i_pis, *i_w, i_target->energy );
    } // end of loop over betas.
    __DEBUGTRYEND(, "Error in Fit::create_A_n_b.\n" )
  } 

  template< class T_POLICY >
  opt::ErrorTuple Fit<T_POLICY> :: check_one( const t_Vector &_ecis,
                                              types::t_unsigned _n, 
                                              bool _verbose )
  {
    __DEBUGTRYBEGIN
    namespace fs = boost::filesystem;
    namespace bblas = boost::numeric::ublas;
    const t_Structures :: value_type &structure = structures[_n];
    const t_Weights :: value_type &weight = weights[_n];
    const t_Pis :: value_type pis = pis[_n];
    const std::string name = fs::path( structure.name ).leaf() ;
    const types::t_real target = structure.energy;
 
    types::t_real predic( bblas::inner_prod( _ecis, pis ) );
    opt::ErrorTuple error( target - predic,  weight );
    if( verbose )
      std::cout << "  structure: " << std::setw(30) << name << "  "
                << "Target: " << std::fixed << std::setw(8) 
                << std::setprecision(2) << target << " "
                << "Separable: " << std::fixed << std::setw(8)
                << std::setprecision(2) << predic << "   "
                << "|Target-Separable| * weight: "
                << std::fixed << std::setw(10) << std::setprecision(3) 
                << error.mean()
                << "\n";
    return error;
    __DEBUGTRYEND(, "Error in Fit::check_one().\n" )
  }

  template< class T_POLICY >
  opt::ErrorTuple Fit<T_POLICY> :: check( const t_Vector &_ecis,
                                bool _training, bool _verbose )
  {
    __DEBUGTRYBEGIN
    opt::ErrorTuple result;
    for(types::t_int n(0); n < structures.size(); ++n, ++i_weight )
    {
      bool found = std::find( excluded.begin(), excluded.end(), n ) != excluded.end();
      if( found xor _training ) continue;

      result += check_one( _reg, _ecis, n )
    }
    return result;
    __DEBUGTRYEND(, "Error in Fit::check_all().\n" )
  }

  opt::ErrorTuple Fit :: check( const t_Vector &_ecis,
                                bool _training, bool _verbose )
  {
    __DEBUGTRYBEGIN
    opt::ErrorTuple result;
    for(types::t_int n(0); n < structures.size(); ++n, ++i_weight )
    {
      bool found = std::find( excluded.begin(), excluded.end(), n ) != excluded.end();
      if( found xor _training ) continue;

      result += check_one( _reg, _ecis, n )
    }
    return result;
    __DEBUGTRYEND(, "Error in Fit::check_all().\n" )
  }

}

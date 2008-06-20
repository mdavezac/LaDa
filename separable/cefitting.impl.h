//
//  Version: $Id$
//
namespace Fitting
{
  template< class T_ALLSQ, class T_COLLAPSE >
  types::t_real SepCeInterface :: fitexcept( T_ALLSQ &_allsq,
                                             T_COLLAPSE &_collapse ) const
  {
    // First creates the input vector.
    std::vector<t_Configurations::value_type::first_type> input;
    typename T_ALLSQ :: t_Vector w, y;
    input.reserve( training.size() * training[0].size() );
    w.reserve( training.size() * training[0].size() );
    y.reserve( training.size() * training[0].size() );
    std::vector< t_Configurations > :: const_iterator i_train( training.begin() );
    std::vector< t_Configurations > :: const_iterator i_train_end( training.end() );
    std::vector< types::t_real > :: const_iterator i_weight( weight.begin() );
    std::vector< types::t_real > :: const_iterator i_target( targets.begin() );
    for(types::t_unsigned i=0; i_train != i_train_end ;
        ++i, ++i_train, ++i_weight, +i_target )
    {
      if(     exclude.size() 
          and std::find( exclude.begin, exclude.end(), i ) == exclude.end() ) 
        continue;
      t_Configurations :: const_iterator i_conf( i_train->begin() );
      t_Configurations :: const_iterator i_conf_end( i_train->end() );
      for(; i_conf != i_conf_end; ++i_conf )
      {
        input.append( i_conf->second ); 
        w.resize( weight.size() + i_conf->second->size(),
                  (*i_weight) * i_conf->first );
        y.resize( y.size() + i_conf->second->size(), 
                  *i_target );
      }
    }
    // initializes the collapse functor.
    _collapse.reset();
    _collapse.init( input );
    allsq.lsq.init_w( _w );
    allsq.lsq.do_weights = true;

    // Then creates the vectors of coefficients with random initial value.
    typename T_ALLSQ :: t_Vectors coefs;
    _collapse.create_coefs( coefs );
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs = coefs.begin();
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs_end = coefs.end();
    for(; i_coefs != i_coefs_end; ++i_coefs )
      std::generate( i_coefs-begin(), i_coefs->end(),
                     boost::ref( rng ) );

    // finally performs fit
    types::t_real convergence = _allsq( *coefs, &_collapse );

    _collapse.reassign( *coefs );
    delete coefs;
  }

  template< class T_ALLSQ, class T_COLLAPSE>
    std::pair< std::pair< types::t_real, types::t_real> >
      leave_one_out( SepCeInterface &_interface,
                    T_ALLSQ &_allsq, 
                    T_COLLAPSE &_collapse, 
                    bool _verbose = false )
    {
      std::pair<types::t_real, types::t_real> training(0,0);
      std::pair<types::t_real, types::t_real> prediction(0,0);
      std::pair<types::t_real, types::t_real> intermediate;
      _interface.exclude.resize(1, 0);
      types::t_unsigned N = _interface.training_set_size();

      for(; _interface.excluded[0] < N; ++_interface.excluded[0] )
      {
        _interface.fit( _allsq, _collapse );

        if( verbose ) std::cout << "Training:\n";
        intermediate = _interface.check_training( _verbose );
        training.first += intermediate.first * (types::t_real) (N-1);
        if( intermediate.second > training.second ) 
          training.second = intermediate.second;

        if( verbose ) std::cout << "Prediction:\n";
        intermediate = _interface.check_training( _verbose );
        prediction.first += intermediate.first;
        if( intermediate.second > prediction.second ) 
          training.second = prediction.second;
      }

      training.first /= (types::t_real) (N*N - N);
      prediction.first /= (types::t_real) N;
      return std::pair< std::pair< types::t_real,
                                   types::t_real > >( training, prediction);
    }
}

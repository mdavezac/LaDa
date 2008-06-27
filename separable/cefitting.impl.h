//
//  Version: $Id$
//

#include<opt/random.h>
#include<algorithm>

namespace Fitting
{
  template< class T_ALLSQ, class T_COLLAPSE >
  void SepCeInterface :: fit( T_ALLSQ &_allsq, T_COLLAPSE &_collapse ) const
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
        ++i, ++i_train, ++i_weight, ++i_target )
    {
      if(     exclude.size() 
          and std::find( exclude.begin(), exclude.end(), i ) == exclude.end() ) 
        continue;

      // all energies of equivalent structures are ... equivalent.
      // hence resize rather than pushback.
      y.resize( y.size() + i_train->size(),  *i_target );
      t_Configurations :: const_iterator i_conf( i_train->begin() );
      t_Configurations :: const_iterator i_conf_end( i_train->end() );
      for(; i_conf != i_conf_end; ++i_conf )
      {
        input.push_back( i_conf->first ); 
        w.push_back( (*i_weight) * i_conf->second );
      }
    }
    // initializes the collapse functor.
    _collapse.reset();
    _collapse.init( input, w );
    _allsq.init_targets( y );
    _allsq.llsq.doweights = false;

    // Then creates the vectors of coefficients with random initial value.
    typename T_ALLSQ :: t_Vectors coefs;
    _collapse.create_coefs( coefs );
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs = coefs.begin();
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs_end = coefs.end();
    for(; i_coefs != i_coefs_end; ++i_coefs )
      std::for_each
      (
        i_coefs->begin(), i_coefs->end(),
        boost::lambda::_1 = 3e0-boost::lambda::bind( &opt::random::rng ) - 1.5e0
      );

    // finally performs fit
    types::t_real convergence = _allsq( coefs, &_collapse );

    _collapse.reassign( coefs );
  }

  template< class T_ALLSQ, class T_COLLAPSE>
    std::pair< SepCeInterface::t_PairErrors, SepCeInterface::t_PairErrors> 
      leave_one_out( SepCeInterface &_interface,
                    T_ALLSQ &_allsq, 
                    T_COLLAPSE &_collapse, 
                    bool _verbose = false )
    {
      typedef SepCeInterface::t_PairErrors t_PairErrors;
      typedef std::pair< t_PairErrors, t_PairErrors > t_Pairs;
      t_PairErrors training(0,0);
      t_PairErrors prediction(0,0);
      t_PairErrors intermediate;
      _interface.exclude.resize(1, 0);
      types::t_unsigned N = _interface.training_set_size();

      for(; _interface.exclude[0] < N; ++_interface.exclude[0] )
      {
        _interface.fit( _allsq, _collapse );

        if( _verbose ) std::cout << "Training:\n";
        intermediate = _interface.check_training( _collapse.function, _verbose );
        training.first += intermediate.first * (types::t_real) (N-1);
        if( intermediate.second > training.second ) 
          training.second = intermediate.second;

        if( _verbose ) std::cout << "Prediction:\n";
        intermediate = _interface.check_training( _collapse.function, _verbose );
        prediction.first += intermediate.first;
        if( intermediate.second > prediction.second ) 
          training.second = prediction.second;
      }

      training.first /= (types::t_real) (N*N - N);
      prediction.first /= (types::t_real) N;
      return t_Pairs( training, prediction);
    }
}

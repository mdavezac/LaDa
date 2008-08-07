//
//  Version: $Id$
//

#include <utility>
#include <algorithm>
#include "bestof.h"

namespace Fitting
{
  template< class T_ALLSQ, class T_COLLAPSE >
  types::t_real SepCeInterface :: fit( T_ALLSQ &_allsq, T_COLLAPSE &_collapse ) const
  {
    // First creates the input vector.
    std::vector<t_Configurations::value_type::first_type> input;
    typename std::vector< typename T_ALLSQ :: t_Vector :: value_type > w, y;
    typename std::vector< std::vector< types::t_real > > eweights;
    input.reserve( training.size() * training[0].size() );
    // w.reserve( training.size() * training[0].size() );
    y.reserve( training.size() * training[0].size() );
    std::vector< t_Configurations > :: const_iterator i_train( training.begin() );
    std::vector< t_Configurations > :: const_iterator i_train_end( training.end() );
    std::vector< types::t_real > :: const_iterator i_weight( weights.begin() );
    std::vector< types::t_real > :: const_iterator i_target( targets.begin() );
    for(types::t_unsigned i=0; i_train != i_train_end ;
        ++i, ++i_train, ++i_weight, ++i_target )
    {
      if(     exclude.size() 
          and std::find( exclude.begin(), exclude.end(), i ) != exclude.end() ) 
        continue;

      y.push_back( *i_target + offset );
      w.push_back( *i_weight );
      t_Configurations :: const_iterator i_conf( i_train->begin() );
      t_Configurations :: const_iterator i_conf_end( i_train->end() );
      std::vector< types::t_real > dummy;
      for(; i_conf != i_conf_end; ++i_conf )
      {
        input.push_back( i_conf->first ); 
        dummy.push_back( i_conf->second );
      }
      eweights.push_back( dummy );
    }

    // initializes the collapse functor.
    _collapse.reset();
    _collapse.init( input, w, eweights );
    _allsq.init_targets( y );
    w.clear();

    // Then creates the vectors of coefficients with random initial value.
    typename T_ALLSQ :: t_Vectors coefs;
    _collapse.create_coefs( coefs, howrandom );
    if( nb_initial_guesses > 1 )
      Separables::bestofcoefs( coefs, y, _collapse, howrandom,
                               nb_initial_guesses, verbose );
    
    // finally performs fit
    types::t_real residual = _allsq( coefs, &_collapse );

    _collapse.reassign( coefs );
    return residual;
  }

  template< class T_ALLSQ, class T_COLLAPSE>
    std::pair< opt::ErrorTuple, opt::ErrorTuple > 
      leave_one_out( SepCeInterface &_interface,
                    T_ALLSQ &_allsq, 
                    T_COLLAPSE &_collapse, 
                    bool _verbose )
    {
      opt::ErrorTuple training, prediction;
      types::t_unsigned N = _interface.training_set_size();

      bool first_iter = true;
      _interface.exclude.resize(1);
      for(; _interface.exclude[0] < N; ++_interface.exclude[0], first_iter=false )
      {
        opt::ErrorTuple intermediate;
        _interface.fit( _allsq, _collapse );

        if( _verbose ) std::cout << "Training:\n";
        intermediate = _interface.check_training( _collapse.function, _verbose );
        if( _verbose ) std::cout << intermediate << "\n";
        training += opt::ErrorTuple( intermediate.mean(), 1e0 );

        if( _verbose ) std::cout << "Prediction:\n";
        intermediate = _interface.check_predictions( _collapse.function, _verbose );
        if( _verbose ) std::cout << intermediate << "\n";
        prediction += intermediate;
      }

      return std::make_pair( training, prediction );
    }
}

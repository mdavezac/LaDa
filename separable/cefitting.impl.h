//
//  Version: $Id$
//
namespace Fitting
{
  template< class T_ALLSQ, class T_COLLAPSE = Separable::Collapse >
  types::t_real SepCeInterface :: fitexcept( T_ALLSQ &_allsq, T_COLLAPSE *_collapse,
                                             std::vector< types::t_unsigned > *_excpt=NULL) const
  {
    // First creates the input vector.
    std::vector<t_Configurations::value_type::first_type> input;
    typename T_ALLSQ :: t_Vector w;
    std::vector< t_Configurations > :: const_iterator i_train( training.begin() );
    std::vector< t_Configurations > :: const_iterator i_train_end( training.end() );
    std::vector< types::t_real > :: const_iterator i_weight = weight.begin();
    for(types::t_unsigned u=0; i_train != i_train_end ; ++i, ++i_train, ++i_weight )
    {
      if(     _excpt
          and std::find( _excpt->begin, _excpt->end(), i ) == _excpt.end() ) 
        continue;
      t_Configurations :: const_iterator i_conf( i_train->begin() );
      t_Configurations :: const_iterator i_conf_end( i_train->end() );
      for(; i_conf != i_conf_end; ++i_conf )
      {
        input.append( i_conf->second ); 
        weight.resize( weight.size() + i_conf->second->size(),
                       (*i_weight) * i_conf->first );
      }
    }
    // initializes the collapse functor.
    collapse->reset();
    collapse->init( input );

    // Then creates the vectors of coefficients with random initial value.
    typename T_ALLSQ :: t_Vectors coefs;
    collapse->create_coefs( coefs );
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs = coefs.begin();
    typename T_ALLSQ :: t_Vectors :: iterator i_coefs_end = coefs.end();
    for(; i_coefs != i_coefs_end; ++i_coefs )
      std::generate( i_coefs-begin(), i_coefs->end(),
                     boost::ref( rng ) );

    // finally performs fit
    types::t_real convergence = _allsq( *coefs, collapse );

    if( not _excpt ) collapse->reassign( *coefs );
    delete coefs;
  }

}

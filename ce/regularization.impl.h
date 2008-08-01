//
//  Version: $Id$
//

namespace CE
{
  template< class T_MINIMIZER >
    void drautz_diaz_ortiz( Regulated &_reg, 
                            const T_MINIMIZER &_minimizer,
                            types::t_int _verbosity,
                            types::t_real _initweights )
    {
      __DEBUGTRYBEGIN
      namespace bl = boost::lambda;
 
      std::cout << "Starting Variational figure plot.\n";
 
      Regulated::t_Clusters save_clusters( _reg.clusters );
 
      
      types::t_unsigned nb_cls( _reg.clusters.size() );
 
      while( _reg.clusters.size() > 0 )
      {
        Regulated :: t_Vector solution( nb_cls );
        Regulated :: t_Vector ecis( nb_cls );
        Regulated :: t_Vector zero_vec( nb_cls, 0e0);
        const types::t_real range(100);
        std::for_each
        ( 
          solution.begin(), solution.end(),
          bl::_1 =  _initweights //bl::bind( &opt::random::rng ) * range - range * 0.5e0
        );
        std::for_each
        ( 
          ecis.begin(), ecis.end(),
          bl::_1 =  bl::bind( &opt::random::rng ) * range - range * 0.5e0
        );
        // Fitting Error
        types::t_real fit = _reg.fit( ecis );
 
        // CV with zero weights
        types::t_real cvwz( _reg( &zero_vec[0] ) ); 
 
        // CV with optimized weights
        types::t_real cvw;
        try { cvw = _minimizer( _reg, solution ); }
        catch( std::exception &_e )
          { std::cout << __SPOT_ERROR << _e.what(); }
        
        // reduces figures by one.
        _reg.reassign( solution );
        if( _verbosity >= 5 )
        {
          Regulated::t_Clusters :: const_iterator i_clusters = _reg.clusters.begin();
          Regulated::t_Clusters :: const_iterator i_clusters_end = _reg.clusters.end();
          for(; i_clusters != i_clusters_end; ++i_clusters )
            std::cout << i_clusters->front();
        }
        Cluster cluster( _reg.reduce() );
        --nb_cls;
 
        std::cout << " Number of clusters: " << nb_cls << "\n"
                  << "  CV with weights: " << cvw << "\n"
                  << "  CV with weights=0: " << cvwz << "\n"
                  << "  Fitting squared error: " << fit << "\n"
                  << "  Dropping " << cluster << "\n\n";
      }
      __DEBUGTRYEND(, "Error encountered in procedure drautz_ortiz_diaz.\n" )
    }
}

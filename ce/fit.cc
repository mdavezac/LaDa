//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iomanip>
#include<algorithm>
#include<numeric>
#include<iostream>
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/filesystem/path.hpp>
#include<boost/filesystem/operations.hpp>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>

#include "fit.h"

namespace LaDa
{
namespace CE
{
  opt::ErrorTuple BaseFit :: check_one( const t_Vector &_ecis,
                                        types::t_unsigned _n, 
                                        bool _verbose ) const
  {
    __DEBUGTRYBEGIN
    namespace fs = boost::filesystem;
    namespace bblas = boost::numeric::ublas;
    const t_Structures :: value_type &structure = structures()[_n];
    const t_Weights :: value_type &weight = weights()[_n];
    const std::string name = fs::path( structure.name ).leaf() ;
    const types::t_real target = structure.energy;
 
    types::t_real predic( bblas::inner_prod( _ecis, pis[_n] ) );
    opt::ErrorTuple error( target - predic,  weight );
    if( _verbose )
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

  void BaseFit :: reassign( const t_Vector &_arg, t_Clusters &_clusters ) const
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
  
    // now copies to clusters.
    t_Clusters :: iterator i_clusters = _clusters.begin();
    t_Clusters :: iterator i_clusters_end = _clusters.end();
    t_Vector :: const_iterator i_arg = _arg.begin();
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_arg )
      std::for_each
      (
        i_clusters->begin(), i_clusters->end(),
        bl::bind( &Cluster::eci, bl::_1 ) = bl::constant(*i_arg)
      );
    __DEBUGTRYEND(, "Error in Regulated::reassign().\n" )
  }

  void BaseFit :: init( const t_Clusters &_clusters )
  {
    find_pis( _clusters, structures(), pis );
    nb_cls = _clusters.size();
    weights().resize( structures().size() );
    std::fill( weights().begin(), weights().end(), 1e0 );
  }

  opt::NErrorTuple BaseFit :: mean_n_var() const
  {
    namespace bl = boost::lambda;
    t_Weights targets( structures().size() );
    std::transform
    (
      structures().begin(), structures().end(), targets.begin(),
      bl::bind( &Crystal::Structure::energy, bl::_1 )
    );
    __ASSERT( weights().size() != targets.size(),
              "Inconsistent weight and target sizes.\n" )
    return opt::mean_n_var( targets, weights() );
  }

} // end of namespace CE
}

//
//  Version: $Id$
//

#ifndef _LADA_CE_DRAUTZ_DIAZ_ORTIZ_H_
#define _LADA_CE_DRAUTZ_DIAZ_ORTIZ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <crystal/structure.h>
#include <minimizer/cgs.h>
#include <minimizer/variant.h>

#include "regularization.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Computes CV scores and reduces number of clusters to zero.
    //! \details Regulated::clusters are unchanged at the end of the run.
    //! \brief Regulated Cluster-Expansion.
    //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
    template< class T_MPLVECTOR >
      void drautz_diaz_ortiz( Regulated &_reg, 
                              const LaDa::Minimizer::Variant<T_MPLVECTOR> &_minimizer,
                              types::t_int _verbosity,
                              types::t_real _initweights )
      {
        __DEBUGTRYBEGIN
        namespace bl = boost::lambda;
   
        std::cout << "Starting Variational figure plot.\n";
   
        Regulated::t_Clusters save_clusters( _reg.clusters );
   
        
        types::t_unsigned nb_cls( _reg.clusters.size() );
   
        while( _reg.clusters.size() > 2 )
        {
          std::cout << " Number of clusters in optimization set: " << nb_cls << "\n";
          if( _verbosity >= 2 )
          {
            std::cout << " List of Clusters:\n";
            Regulated::t_Clusters :: const_iterator i_clusters = _reg.clusters.begin();
            Regulated::t_Clusters :: const_iterator i_clusters_end = _reg.clusters.end();
            for(; i_clusters != i_clusters_end; ++i_clusters )
              std::cout << i_clusters->front();
          }
          Regulated :: t_Arg solution( nb_cls );
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
          opt::ErrorTuple fitnoreg = _reg.fit( ecis, &zero_vec[0] );
   
          // CV with zero weights
          types::t_real cvwz( _reg( zero_vec ) ); 
   
          // CV with optimized weights
          types::t_real cvw;
          __TRYBEGIN
           cvw = _minimizer( _reg, solution ); 
          __TRYEND(," ") 
          std::for_each( solution.begin(), solution.end(), std::cout << bl::_1 << " " );
          std::cout << "\n";
          
          // Fitting Error
          opt::ErrorTuple fitwithreg = _reg.fit( ecis, &solution[0] );
          std::for_each( ecis.begin(), ecis.end(), std::cout << bl::_1 << " " );
          std::cout << "\n";
          
          // reduces figures by one.
          _reg.reassign( ecis );
          Cluster cluster( _reg.reduce() );
          --nb_cls;
   
          std::cout << "  CV with weights: " << cvw << "\n"
                    << "  CV with weights=0: " << cvwz << "\n"
                    << "  Fitting squared error (no regulation): " << fitnoreg << "\n"
                    << "  Fitting squared error (with regulation): " << fitwithreg << "\n"
                    << "  Dropping " << cluster << "\n\n";
        }
        __DEBUGTRYEND(, "Error encountered in procedure drautz_ortiz_diaz.\n" )
      }
  }
} // namespace LaDa

#endif

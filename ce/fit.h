//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifndef _CE_FIT_H_
#define _CE_FIT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/cgs.h>
#include <opt/gsl_simplex.h>
#include <opt/errors.h>
#include <crystal/structure.h>

#include "cluster.h"

namespace CE
{
  //! A class for fitting cluster-expansions to LDA data.
  template< class T_POLICY >
  class Fit : public T_POLICY
  {
    public:
      //! Type of the fitting matrices.
      typedef boost::numeric::ublas::matrix<types::t_real> t_Matrix;
      //! Type of the fitting target vectors.
      typedef boost::numeric::ublas::vector<types::t_real> t_Vector;
      //! A container of structures.
      typedef std::vector< Crystal::Structure > t_Structures;

    protected:
      //! Type of the class representing a cluster.
      typedef Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef std::vector< t_Cluster > t_EquivClusters;
      //! A container of Pis for a single structure.
      typedef t_Vector t_StructPis;
      //! A container of weights.
      typedef std::vector< types::t_real > t_Weights;

    public:
      //! Container of classes of equivalent clusters.
      typedef std::vector< t_EquivClusters > t_Clusters;
      //! A container of Pis for a single structure.
      typedef std::vector< t_StructPis > t_Pis;

    public:
      //! A container of structures.
      t_Structures structures;
      
      //! Constructor.
      Fit() {};
      //! Destructor.
      ~Fit() {};

      //! Evaluates the cv score for the weights on input.
      opt::ErrorTuple operator()( t_Vector &_arg ) const; 

      //! Initializes from structures. 
      void init( const t_Clusters &_clusters )
      {
        find_pis( _clusters, structures, pis );
        nb_cls = _clusters.size();
        weights.resize( structures.size() );
        std::fill( weights.begin(), weights.end(), 1e0 );
      }

      //! Reassigns ecis from argument values.
      void reassign( const t_Vector &_arg ) const;

      //! Computes error for training set.
      opt::ErrorTuple check_training( const t_Vector &_ecis,
                                      bool _verbose ) const
       { return check( _ecis, true, _verbose ); }
      //! Computes error for prediction set.
      opt::ErrorTuple check_prediction( const t_Vector &_ecis,
                                      bool _verbose ) const
       { return check( _ecis, true, _verbose ); }

    protected:
      //! Computes \a _A and \a _b excluding excluded structures.
      void create_A_n_b( t_Vector &_A, t_Vector &_b ) const;
      //! Adds current pis to \a _A and \a _b.
      void add_to_A_n_b( t_Vector &_A, t_Vector &_b,
                         const t_StructPis &_pis,
                         const types::t_real _weight,
                         const types::t_real _energy ) const;
      //! Computes error for training or prediction set.
      opt::ErrorTuple check( const t_Vector &_ecis,
                             bool _training,
                             bool _verbose ) const
      //! Computes error for one structure.
      opt::ErrorTuple check_one( const t_Vector &_ecis,
                                 types::t_unsigned _n, 
                                 bool _verbose = false ) const;

      //! A container of pis for all structures.
      t_Pis pis;
      //! A container of weights.
      t_Weights weights;
  };

  //! A class for fitting cluster-expansions to LDA data, excluding some.
  class ExFit : public Fit
  {
    public:
      //! Type of the fitting matrices.
      typedef Fit :: t_Matrix t_Matrix;
      //! Type of the fitting target vectors.
      typedef Fit :: t_Vector t_Vector;
      //! A container of structures.
      typedef Fit :: t_Structures t_Structures;

    protected:
      //! Type of the class representing a cluster.
      typedef Fit :: t_Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef Fit :: t_EquivClusters t_EquivClusters;
      //! A container of Pis for a single structure.
      typedef Fit :: t_StructPis t_StructPis;
      //! A container of weights.
      typedef Fit :: t_Weights t_Weights;

    public:
      //! Container of classes of equivalent clusters.
      typedef Fit :: t_Clusters t_Clusters;
      //! A container of Pis for a single structure.
      typedef Fit :: t_Pis t_Pis;
      //! A container of Pis for a single structure.
      typedef std::vector< types::t_unsigned > t_Excluded;

    public:
      //! A container of structures.
      t_Structures structures;
      //! \brief A container of indices to the structures which should be
      //!        excluded from the fit.
      t_Excluded excluded;
      
      //! Constructor.
      Fit() {};
      //! Destructor.
      ~Fit() {};

      //! Evaluates the cv score for the weights on input.
      opt::ErrorTuple operator()( t_Vector &_arg ) const; 

      //! Initializes from structures. 
      void init( const t_Clusters &_clusters )
      {
        find_pis( _clusters, structures, pis );
        nb_cls = _clusters.size();
        weights.resize( structures.size() );
        std::fill( weights.begin(), weights.end(), 1e0 );
      }

      //! Reassigns ecis from argument values.
      void reassign( const t_Vector &_arg ) const;

      //! Computes error for training set.
      opt::ErrorTuple check_training( const t_Vector &_ecis,
                                      bool _verbose ) const
       { return check( _ecis, true, _verbose ); }
      //! Computes error for prediction set.
      opt::ErrorTuple check_training( const t_Vector &_ecis,
                                      bool _verbose ) const
       { return check( _ecis, true, _verbose ); }

    protected:
      //! Computes \a _A and \a _b excluding excluded structures.
      void create_A_n_b( t_Vector &_A, t_Vector &_b ) const;
      //! Adds current pis to \a _A and \a _b.
      void add_to_A_n_b( t_Vector &_A, t_Vector &_b,
                         const t_StructPis &_pis,
                         const types::t_real _weight,
                         const types::t_real _energy ) const;
      //! Computes error for training or prediction set.
      opt::ErrorTuple Fit :: check( const t_Vector &_ecis,
                                    bool _training,
                                    bool _verbose ) const
      //! Computes error for one structure.
      opt::ErrorTuple check_one( const t_Vector &_ecis,
                                 types::t_unsigned _n, 
                                 bool _verbose = false ) const;

      //! A container of pis for all structures.
      t_Pis pis;
      //! A container of weights.
      t_Weights weights;
  };

} // end of namespace CE

#include "regularization.impl.h"

#endif 

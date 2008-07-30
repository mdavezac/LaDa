//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifndef _CE_REGULARIZATION_H_
#define _CE_REGULARIZATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/cgs.h>
#include <crystal/structure.h>

#include "cluster.h"

namespace CE
{
  // forward declaration
  //! \cond
  class Regulated;
  //! \endcond

  //! \brief Computes CV scores and reduces number of clusters to zero.
  //! \details Regulated::clusters are unchanged at the end of the run.
  //! \brief Regulated Cluster-Expansion.
  //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
  void drautz_diaz_ortiz( Regulated &_reg,
                          types::t_real tolerance = 1e-4,
                          types::t_int _verbosity = 0 );

  //! \brief Regulated Cluster-Expansion.
  //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
  //!      and Alejandro Diaz-Ortiz, PRB \bf 73, 224207 (2007)</A>.
  class Regulated
  {
    mutable bool doinggradient;
    public:
      //! Type of the fitting matrices.
      typedef boost::numeric::ublas::matrix<types::t_real> t_Matrix;
      //! Type of the fitting target vectors.
      typedef boost::numeric::ublas::vector<types::t_real> t_Vector;
      //! A container of structures.
      typedef std::vector< Crystal::Structure > t_Structures;
      //! Type of the return.
      typedef types::t_real t_Return;
      //! Type of the input variables.
      typedef boost::numeric::ublas::vector<t_Return> t_Arg;

    protected:
      //! Type of a pair of fitting matrix and vector.
      typedef std::pair< t_Matrix, t_Vector > t_FittingPair;
      //! A container of fitting pairs.
      typedef std::vector< t_FittingPair > t_FittingPairs;
      //! Type of the ecis.
      typedef std::vector< t_Vector > t_FittedEcis;
      //! Type of the class representing a cluster.
      typedef Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef std::vector< t_Cluster > t_EquivClusters;
      //! A container of Pis for a single structure.
      typedef t_Vector t_StructPis;
      //! A container of weights.
      typedef std::vector< types::t_real > t_Weights;
      //! A container of targets.
      typedef std::vector< types::t_real > t_Targets;

    public:
      //! Container of classes of equivalent clusters.
      typedef std::vector< t_EquivClusters > t_Clusters;
      //! A container of Pis for a single structure.
      typedef std::vector< t_StructPis > t_Pis;

    public:
      //! The clusters to fit.
      t_Clusters clusters;
      //! The fitting procedure.
      Fitting::Cgs cgs;


      //! Constructor.
      Regulated() {};
      //! Destructor.
      ~Regulated() {};

      //! Evaluates the cv score for the weights on input.
      t_Return operator()( const types::t_real * _arg ) const; 
      //! Evaluates the gradient.
      void gradient( const types::t_real * _arg,
                     types::t_real *_gradient ) const;
      //! Initializes from structures. 
      void init( const t_Structures &_structures );
      //! Reduce regulated function and cluster by 1.
      types::t_unsigned reduce();
      //! Fit to all data using linear-least square fit.
      types::t_real fit( t_Vector &_x );
      //! Analytical fit
      types::t_real anafit();
      //! Computes square errors.
      types::t_unsigned square_errors() const;
      //! Reassigns ecis from argument values.
      void reassign( const t_Arg &_arg );

    protected:
      //! Initializes Regulated::sums.
      void init_sums();
      //! \brief Constructs \a _A and \a _b fitting matrix and vector excluding
      //!        structure \a _k.
      //! \details Does not include wights.
      void construct_pair( t_FittingPair &_pair, types::t_unsigned &_k );
      //! Constructs all fitting pairs, without the weights.
      void construct_pairs();

    protected:
      //! The number of clusters.
      types::t_unsigned nb_cls;
      //! A container of pis for all structures.
      t_Pis pis;
      //! Type of Regulated::esums.
      typedef t_Vector t_ESums;
      //! Type of Regulated::psums.
      typedef t_Matrix t_PSums;
      //! \f$=\sum_s w_s \phi_{\alpha, s}\phi_{\beta, s}\f$.
      t_PSums psums;
      //! \f$=\sum_s w_s E_s\phi_{\beta, s}\f$.
      t_ESums esums;
      //! A container of weights.
      t_Weights weights;
      //! A container of weights.
      t_Targets targets;
      //! A container of fitting matrices and vectors.
      t_FittingPairs fittingpairs;
      //! A container of fitted interactions.
      mutable t_FittedEcis fittedecis;
  };

} // end of namespace CE

#endif 

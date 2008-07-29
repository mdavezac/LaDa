//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifndef _CE_REGULARIZATION_H_
#define _CE_REGULARIZATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
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
    public:
      //! Type of the class representing a cluster.
      typedef Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef std::vector< t_Cluster > t_EquivClusters;
      //! Container of classes of equivalent clusters.
      typedef std::vector< t_EquivClusters > t_Clusters;
      //! A container of Pis for a single structure.
      typedef std::vector< types::t_real > t_StructPis;
      //! A container of Pis for a single structure.
      typedef std::vector< t_StructPis > t_Pis;
      //! A container of structures.
      typedef std::vector< Crystal::Structure > t_Structures;
      //! A container of weights.
      typedef std::vector< types::t_real > t_Weights;
      //! A container of targets.
      typedef std::vector< types::t_real > t_Targets;
      //! Type of the Jacobian.
      typedef std::vector<types::t_real> t_Jacobian;
      //! Type of the return.
      typedef std::vector<types::t_real> t_Return;
      //! Type of the input variables.
      typedef std::vector< types::t_real > t_Arg;


    public:
      //! The clusters to fit.
      t_Clusters clusters;

      //! Constructor.
      Regulated() {};
      //! Destructor.
      ~Regulated() {};

      //! Evaluates the cv score for the weights on input.
      void operator()( const types::t_real * _arg,
                       types::t_real *_return ) const; 
      //! Evaluates the jacobian.
      void jacobian( const types::t_real * _arg,
                     types::t_real *_jacobian ) const;
      //! Initializes from structures. 
      void init( const t_Structures &_structures );
      //! Dimension of the function.
      types::t_unsigned dim() const { return targets.size(); }
      //! Reduce regulated function and cluster by 1.
      types::t_unsigned reduce();
      //! Fit to all data using linear-least square fit.
      types::t_real fit( types::t_real _tol = 1e-4, 
                         types::t_unsigned _imax = 40, 
                         bool verbose = false );
      //! Computes square errors.
      types::t_unsigned square_errors() const;
      //! Reassigns ecis from argument values.
      void reassign( const t_Arg &_arg );

    protected:
      //! Initializes Regulated::sums.
      void init_sums();

    protected:
      //! The number of clusters.
      types::t_unsigned nb_cls;
      //! A container of pis for all structures.
      t_Pis pis;
      //! Type of Regulated::numsums.
      typedef std::vector< types::t_real > t_NumSums;
      //! Type of Regulated::denumsums.
      typedef std::vector< std::vector< types::t_real > > t_DenumSums;
      //! \f$=\sum_s w_s \phi_{\alpha, s}\phi_{\beta, s}\f$.
      t_DenumSums denumsums;
      //! \f$=\sum_s w_s E_s\phi_{\beta, s}\f$.
      t_NumSums numsums;
      //! A container of weights.
      t_Weights weights;
      //! A container of weights.
      t_Weights targets;
  };

} // end of namespace CE

#endif 

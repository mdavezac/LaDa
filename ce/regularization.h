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
#include <boost/tuple/tuple.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/cgs.h>
#include <opt/gsl_simplex.h>
#include <crystal/structure.h>

#include "cluster.h"

namespace CE
{
  // forward declaration
  //! \cond
  class Regulated;
  class Fit;
  namespace details
  {
    //! An error tuple.
    class ErrorTuple : public boost::tuple<types::t_real, types::t_real, types::t_real> 
    {  
      //! Base class.
      typedef boost::tuple<types::t_real, types::t_real, types::t_real> t_Base;
      public:
        ErrorTuple() { get<0>() = 0e0; get<1>() = 0e0; get<2>() = 0e0; }
        ErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                  : t_Base( _a, _b, _c ) {} 
        ErrorTuple  ( types::t_real _a, types::t_real _b )
                  : t_Base( _a * _a * _b, _a * _b, std::max( _a, get<2>() ) ) {} 
        ErrorTuple  ( types::t_real _a )
                  : t_Base( _a * _a, _a, std::max( _a, get<2>() ) ) {} 
        ErrorTuple  ( const ErrorTuple &_e ) : t_Base( _e ) {} 
    };
    //! A normalized error tuple.
    struct NErrorTuple : public ErrorTuple
    {
      //! Base class.
      typedef ErrorTuple t_Base;
      public:
        //! Variance.
        types::t_real variance;
        //! Mean.
        types::t_real mean;

        NErrorTuple() : t_Base(), variance(0), mean(0) {}
        NErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                   : t_Base( _a, _b, _c ), variance(0), mean(0) {} 
        NErrorTuple  ( types::t_real _a, types::t_real _b )
                   : t_Base( _a * _a * _b, _a * _b, std::max( _a, get<2>() ) ),
                     variance(0), mean(0) {}
        NErrorTuple  ( types::t_real _a )
                   : t_Base( _a * _a, _a, std::max( _a, get<2>() ) ),
                     variance(0), mean(0) {}
        NErrorTuple  ( const NErrorTuple &_e )
                   : t_Base( _e ), variance( _e.variance ), mean( _e.mean ) {}
        const NErrorTuple& operator=( const ErrorTuple &_e )
        { 
          get<0>() = _e.get<0>(); get<1>() = _e.get<1>(); get<2>() = _e.get<2>(); 
          return *this;
        }
    };
  }
  //! \endcond

  //! \brief Computes CV scores and reduces number of clusters to zero.
  //! \details Regulated::clusters are unchanged at the end of the run.
  //! \brief Regulated Cluster-Expansion.
  //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
  template< class T_MINIMIZER >
  void drautz_diaz_ortiz( Regulated &_reg,
                          const T_MINIMIZER &_minimizer,
                          types::t_int _verbosity = 0,
                          types::t_real _initweights = 0e0 );


  //! \brief Regulated Cluster-Expansion.
  //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
  //!      and Alejandro Diaz-Ortiz, PRB \bf 73, 224207 (2007)</A>.
  class Regulated
  {
    friend class Fit;
    friend void leave_one_out( const Regulated &_reg,
                     const boost::numeric::ublas::vector<types::t_real> &_weights,
                     const std::vector< Crystal::Structure > &_strs );
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


  class Fit 
  {
    public:
      //! Index of structures excluded from the set.
      typedef std::vector< Crystal::Structure > t_Structures;
      //! A container of structures.
      t_Structures structures;
      //! lambda for pair regulation
      types::t_real lambda;
      //! t for pair regulation
      types::t_real tcoef;
      //! Wether to perform pair regulation.
      bool do_pairreg;
      //! Which pair regulation to perform: Laks or Volkers?
      bool laksreg;
      //! Wether to be verbose.
      bool verbose;

      
      //! Constructor.
      Fit() : lambda(0), tcoef(0), do_pairreg(false), laksreg(false), verbose(false) {}
      //! Destructor.
      ~Fit() {};

      //! Performs leave-one-out.
      void leave_one_out( const Regulated &_reg );
      //! Performs leave-one-out.
      void fit( const Regulated &_reg );

    protected:
      //! An error tuple.
      typedef details::ErrorTuple t_ErrorTuple;
      //! Compute pair regulation terms.
      //! Compute pair regulation terms.
      void pair_reg( const Regulated &_reg, Regulated::t_Vector &_weights );
      //! Performs fit (e.g. leave-none-out).
      void fit_but_one( const Regulated &_reg,
                        Regulated :: t_Vector &_x,
                        const Regulated :: t_Vector &_weights,
                        const types::t_unsigned _n ) const;
      //! Check results.
      types::t_real check_one( const Regulated &_reg, 
                               const Regulated :: t_Vector &_ecis,
                               types::t_unsigned _n );
      //! Check all (but one )
      t_ErrorTuple check_all( const Regulated &_reg, 
                              const Regulated :: t_Vector &_ecis,
                              types::t_int _n = -1 );
      //! computes mean and variance of data
      details :: NErrorTuple mean_n_var( const Regulated &_reg ); 
  };

} // end of namespace CE

#include "regularization.impl.h"

#endif 

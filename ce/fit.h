//
//  Version: $Id$
//

#ifndef _CE_FIT_H_
#define _CE_FIT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector> 
#include <utility> 

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <crystal/structure.h>

#include "cluster.h"

namespace CE
{
  // Forward declarations.
  //! \cond 
  class BaseFit;
  template< class T_POLICY > class Fit;
  template< class T_POLICY > class Fit;
  namespace FittingPolicy
  {
    template< class T_FIT > class Nothing;
    template< class T_FIT > class Excluded;
    template< class T_POLICY > class PairReg;
  }
  //! \endcond
  
  //! A class for fitting cluster-expansions to LDA data.
  class BaseFit 
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
      //! A container of weights.
      t_Weights weights;
      
      //! Constructor.
      BaseFit() : nb_cls(0) {};
      //! Destructor.
      ~BaseFit() {};

      //! Initializes from structures. 
      void init( const t_Clusters &_clusters );
      //! Reassigns ecis to clusters.
      void reassign( const t_Vector &_arg, t_Clusters &_clusters ) const;
      //! Computes mean and variance.
      opt::NErrorTuple mean_n_var() const;

    protected:
      //! Computes error for one structure.
      opt::ErrorTuple check_one( const t_Vector &_ecis,
                                 types::t_unsigned _n, 
                                 bool _verbose = false ) const;

      //! A container of pis for all structures.
      t_Pis pis;
      //! Number of clusters.
      types::t_unsigned nb_cls;
  };

  //! \cond
  namespace details
  {
    template< class T_SOLVER, class T_CLASS >
      opt::ErrorTuple operator_( const T_CLASS &_class,
                                 BaseFit::t_Vector &_x, 
                                 const T_SOLVER &_solver );
  }
  //! \endcond

  //! Performs leave-one-out. Watch out for request \a _fit interface.
  template< class T_FIT, class T_SOLVER >
    std::pair< opt::ErrorTuple, opt::ErrorTuple >
      leave_one_out( const T_FIT &_fit, 
                     const T_SOLVER &_solver,
                     BaseFit::t_Vector &_x,
                     bool _verbose );
  //! Performs leave-one-out with regularization. Watch out for request \a _fit interface.
  template< class T_FIT, class T_SOLVER >
    std::pair< opt::ErrorTuple, opt::ErrorTuple >
      leave_one_out( const T_FIT &_fit,
                     const T_SOLVER &_solver,
                     BaseFit::t_Vector &_x, 
                     const types::t_real *_w, 
                     bool _verbose );

  //! A class for fitting cluster-expansions to LDA data.
  template< class T_POLICY = FittingPolicy::Nothing<BaseFit> >
  class Fit : public T_POLICY
  {
    template< class T_SOLVER, class T_CLASS > friend 
      opt::ErrorTuple details::operator_( const T_CLASS &_class,
                                          BaseFit::t_Vector &_x, 
                                          const T_SOLVER &_solver );
    public:
      //! Some extra policies.
      typedef T_POLICY t_Policy;
      //! Type of the fitting matrices.
      typedef BaseFit::t_Matrix t_Matrix;
      //! Type of the fitting target vectors.
      typedef BaseFit::t_Vector t_Vector;
      //! A container of structures.
      typedef BaseFit::t_Structures t_Structures;

    protected:
      //! Type of the class representing a cluster.
      typedef BaseFit::t_Cluster t_Cluster;
      //! Container of equivalent clusters.
      typedef BaseFit::t_EquivClusters t_EquivClusters;
      //! A container of Pis for a single structure.
      typedef BaseFit::t_StructPis t_StructPis;
      //! A container of weights.
      typedef BaseFit::t_Weights t_Weights;

    public:
      //! Container of classes of equivalent clusters.
      typedef BaseFit::t_Clusters t_Clusters;
      //! A container of Pis for a single structure.
      typedef BaseFit::t_Pis t_Pis;

    public:
      //! Wether to be verbose.
      bool verbose;

      //! Constructor.
      Fit() : verbose(false) {};
      //! Destructor.
      ~Fit() {};

      //! Evaluates the cv score for the weights on input.
      template< class T_SOLVER>
        opt::ErrorTuple operator()( t_Vector &_x,
                                    const T_SOLVER &_solver ) const
          { return details::operator_( *this, _x, _solver ); }
      //! Initializes from structures. 
      void init( const t_Clusters &_clusters )
        { t_Policy :: init( _clusters ); }

    protected:
      //! Computes \a _A and \a _b excluding excluded structures.
      void create_A_n_b( t_Matrix &_A, t_Vector &_b ) const;
      //! Adds current pis to \a _A and \a _b.
      void add_to_A_n_b( t_Matrix &_A, t_Vector &_b,
                         const t_StructPis &_pis,
                         const types::t_real _weight,
                         const types::t_real _energy ) const;

    protected:  
      using t_Policy :: pis;
      using t_Policy :: weights;
      using t_Policy :: nb_cls;
    public:
      using t_Policy :: structures;
  };

  template< class T_POLICY >
    class RegulatedFit : public Fit<T_POLICY> 
    {
      template< class T_SOLVER, class T_CLASS > friend 
        opt::ErrorTuple details::operator_( const T_CLASS &_class,
                                            BaseFit::t_Vector &_x, 
                                            const T_SOLVER &_solver );
      public:
        //! A policy base class.
        typedef T_POLICY t_Policy;
        //! A base class.
        typedef Fit<t_Policy> t_Base;
        
        //! Constructor.
        RegulatedFit () {};
        //! Destructor..
        ~RegulatedFit () {};

        //! Fit with regularization weights.
        template< class T_SOLVER >
        opt::ErrorTuple operator()( BaseFit::t_Vector &_x,
                                    const types::t_real *_weights,
                                    T_SOLVER &_solver ) const;

      protected:
        //! Adds weight regulation to \a _A.
        void other_A_n_b( BaseFit::t_Matrix &_A,
                          BaseFit::t_Vector &_b ) const;

        //! A pointer to the the weights
        mutable const types::t_real *regweights;

      protected:  
        using t_Base :: pis;
        using t_Base :: weights;
        using t_Base :: nb_cls;
      public:
        using t_Base :: structures;
    };

  namespace FittingPolicy
  {
    template< class T_BASE = BaseFit > 
      class Nothing : public T_BASE
      {
        public:
          //! A fitting class.
          typedef T_BASE t_Base;
          //! Constructor.
          Nothing() {};
          //! Destructor.
          ~Nothing() {};

          //! checks training.
          opt::ErrorTuple check_training( const BaseFit::t_Vector &_ecis,
                                          bool _verbose );

        protected:
          //! Computes \a _A and \a _b excluding excluded structures.
          void other_A_n_b( BaseFit::t_Matrix &_A, BaseFit::t_Vector &_b ) const {};
          //! Returns true if index \a _i is an excluded structures.
          bool found( types::t_unsigned _i ) const { return false; }

        protected:  
          using t_Base :: pis;
          using t_Base :: weights;
          using t_Base :: nb_cls;
        public:
          using t_Base :: structures;
      };

    template< class T_BASE = Nothing<> > 
      class Excluded : public T_BASE
      {
        public:
          //! A fitting class.
          typedef T_BASE t_Base;
          //! A container of indices to excluded structues.
          typedef std::vector< types::t_unsigned > t_Excluded;
          //! A container of indices to excluded structues.
          mutable t_Excluded excluded;

          //! Constructor.
          Excluded() {};
          //! Destructor.
          ~Excluded() {};

          //! Computes error for training set.
          opt::ErrorTuple check_training( const BaseFit::t_Vector &_ecis,
                                          bool _verbose ) const
           { return check( _ecis, true, _verbose ); }
          //! Computes error for prediction set.
          opt::ErrorTuple check_predictions( const BaseFit::t_Vector &_ecis,
                                             bool _verbose ) const
           { return check( _ecis, false, _verbose ); }

        protected:
          //! Computes error for training or prediction set.
          opt::ErrorTuple check( const BaseFit::t_Vector &_ecis,
                                 bool _training, bool _verbose ) const;
          //! Returns true if index \a _i is an excluded structures.
          bool found( types::t_unsigned _i ) const;

        protected:  
          using t_Base :: pis;
          using t_Base :: weights;
          using t_Base :: nb_cls;
        public:
          using t_Base :: structures;
      };

    template< class T_BASE = Excluded<> >
      class PairReg : public T_BASE
      {
        public:
          //! A policy base class class.
          typedef T_BASE t_Base;
          //! lambda for pair regulation
          types::t_real lambda;
          //! t for pair regulation
          types::t_real tcoef;
          //! Wether to perform pair regulation.
          bool do_pairreg;
          //! Which pair regulation to perform: Laks or Volkers?
          bool laksreg;

          //! Constructor.
          PairReg () {};
          //! Destructor..
          ~PairReg () {};

        protected:
          //! Possible initialization stuff.
          void init( const BaseFit::t_Clusters &_clusters );
          //! Adds pair regulation to \a _A.
          void other_A_n_b( BaseFit::t_Matrix &_A,
                            BaseFit::t_Vector &_b ) const;

          //! A container to the the pair terms.
          typedef std::vector< std::pair< size_t, types::t_real > > t_PairWeights;
          //! A container to the the pair terms.
          t_PairWeights pairweights;

        protected:  
          using t_Base :: pis;
          using t_Base :: weights;
          using t_Base :: nb_cls;
        public:
          using t_Base :: structures;
      };

  } // end of namespace FittingPolicy

} // end of namespace CE

#include "fit.impl.h"

#endif 

#ifndef _CE_FIT_H_
#define _CE_FIT_H_

#include "LaDaConfig.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <vector> 
#include <utility> 

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <opt/leave_many_out.h>
#include <crystal/structure.h>

#include "cluster.h"

namespace LaDa
{
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
    template< class T_TRAITS > class MixedApproach;
    template< class T_TRAITS > class CEasCollapse;
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
        //! Constructor.
        BaseFit() : nb_cls(0) 
         { structures_.reset( new t_Structures ); weights_.reset( new t_Weights ); };
        //! Copy Constructor.
        BaseFit   ( const BaseFit &_c )
                : nb_cls(_c.nb_cls ), pis( _c.pis ), 
                  structures_( _c.structures_ ), weights_( _c.weights_ ) {}
        //! Destructor.
        ~BaseFit() {}

        //! Initializes from structures. 
        void init( const t_Clusters &_clusters );
        //! Reassigns ecis to clusters.
        void reassign( const t_Vector &_arg, t_Clusters &_clusters ) const;
        //! Computes mean and variance.
        opt::NErrorTuple mean_n_var() const;
        //! Returns a reference to the structures.
        t_Structures& structures()
        { 
          LADA_NASSERT( not structures_.get(), "Empty smart pointer.\n" ) 
          return *structures_;
        }
        //! Returns a constant reference to the structures.
        const t_Structures& structures() const
        { 
          LADA_NASSERT( not structures_.get(), "Empty smart pointer.\n" ) 
          return *structures_;
        }
        //! Returns a reference to the weights.
        t_Weights& weights() 
        { 
          LADA_NASSERT( not weights_.get(), "Empty smart pointer.\n" ) 
          return *weights_;
        }
        //! Returns a constant reference to the weights.
        const t_Weights& weights() const 
        { 
          LADA_NASSERT( not weights_.get(), "Empty smart pointer.\n" ) 
          return *weights_;
        }

      protected:
        //! Computes error for one structure.
        opt::ErrorTuple check_one( const t_Vector &_ecis,
                                   types::t_unsigned _n, 
                                   bool _verbose = false ) const;

        //! A container of pis for all structures.
        t_Pis pis;
        //! Number of clusters.
        types::t_unsigned nb_cls;
        //! A container of structures.
        boost::shared_ptr<t_Structures> structures_;
        //! A container of weights.
        boost::shared_ptr<t_Weights> weights_;
        
    };

    //! \cond
    namespace details
    {
      template< class T_SOLVER, class T_CLASS >
        opt::ErrorTuple operator_( const T_CLASS &_class,
                                   BaseFit::t_Vector &_x, 
                                   const T_SOLVER &_solver );
      template< class T_CLASS, class T_MATRIX, class T_VECTOR >
        void construct_( const T_CLASS &_class, T_MATRIX &_A, T_VECTOR &_b );
    }
    //! \endcond

    //! Performs leave-one-out. Watch out for requested \a _fit interface.
    template< class T_FIT, class T_SOLVER >
      std::pair< opt::ErrorTuple, opt::ErrorTuple >
        leave_one_out( const T_FIT &_fit, 
                       const T_SOLVER &_solver,
                       BaseFit::t_Vector &_x,
                       bool _verbose );
    //! Performs leave-one-out with regularization. Watch out for requested \a _fit interface.
    template< class T_FIT, class T_SOLVER >
      std::pair< opt::ErrorTuple, opt::ErrorTuple >
        leave_one_out( const T_FIT &_fit,
                       const T_SOLVER &_solver,
                       BaseFit::t_Vector &_x, 
                       const types::t_real *_w, 
                       bool _verbose );
    //! Performs leave-many-out for CE fits.
    template< class T_FIT, class T_SOLVER >
      opt::t_ErrorPair leave_many_out( Fitting::LeaveManyOut &_lmo,
                                       const T_FIT &_fit,
                                       const T_SOLVER &_solver );

    //! A class for fitting cluster-expansions to LDA data.
    template< class T_POLICY = FittingPolicy::Nothing<BaseFit> >
    class Fit : public T_POLICY
    {
      template< class T_SOLVER, class T_CLASS > friend 
        opt::ErrorTuple details::operator_( const T_CLASS &_class,
                                            BaseFit::t_Vector &_x, 
                                            const T_SOLVER &_solver );

      template< class T_CLASS, class T_MATRIX, class T_VECTOR > friend
        void details::construct_( const T_CLASS &_class, T_MATRIX &_A, T_VECTOR &_b ) ;
      template< class T_TRAITS > friend class MixedApproach;
      template< class T_TRAITS > friend class CEasCollapse;
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
        //! Copy Constructor.
        template< class TT_POLICY >
          Fit   ( const Fit< TT_POLICY > &_c )
              : BaseFit(_c), verbose(_c.verbose) {}
        //! Destructor.
        ~Fit() {}

        //! Evaluates the cv score for the weights on input.
        template< class T_SOLVER>
          opt::ErrorTuple operator()( t_Vector &_x,
                                      const T_SOLVER &_solver ) const
            { return details::operator_( *this, _x, _solver ); }
        //! \brief Creates fitting matrix and target vector.
        //! \details Does not fit. Adds values to pre-existing matrix and vector.
        template< class T_MATRIX, class T_VECTOR >
          void construct( T_MATRIX &_A, T_VECTOR &_b ) const
            { details::construct_( *this, _A, _b ); }
        //! Initializes from structures. 
        void init( const t_Clusters &_clusters )
          { t_Policy :: init( _clusters ); }

        //! Number of degrees of liberty.
        size_t dof() const { return nb_cls; }

      protected:
        //! Computes \a _A and \a _b excluding excluded structures.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const;

      protected:  
        using t_Policy :: pis;
        using t_Policy :: nb_cls;
        using t_Policy :: weights;

      public:
        using t_Policy :: structures;
        //! Allows rebinding of the separables function.
        template< class TT_POLICY > struct rebind
        {
          //! new separable type.
          typedef Fit< TT_POLICY > type;
        };
    };

    template< class T_POLICY >
      class RegulatedFit : public Fit<T_POLICY> 
      {
        template< class T_SOLVER, class T_CLASS > friend 
          opt::ErrorTuple details::operator_( const T_CLASS &_class,
                                              BaseFit::t_Vector &_x, 
                                              const T_SOLVER &_solver );
        template< class T_CLASS, class T_MATRIX, class T_VECTOR > friend
          void details::construct_( const T_CLASS &_class, T_MATRIX &_A, T_VECTOR &_b );
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
          //! \brief Creates fitting matrix and target vector.
          //! \details Does not fit. Adds values to pre-existing matrix and vector.
          template< class T_MATRIX, class T_VECTOR >
            void construct( T_MATRIX &_A, T_VECTOR &_b ) const
              { details::construct_( *this, _A, _b ); }

          //! Adds weight regulation to \a _A.
          template< class T_MATRIX, class T_VECTOR >
            void other_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const;
    
          //! A pointer to the the weights
          mutable const types::t_real *regweights;

        protected:  
          using t_Base :: pis;
          using t_Base :: weights;
          using t_Base :: nb_cls;
        public:
          using t_Base :: structures;

         //! Allows rebinding of the separables function.
         template< class TT_POLICY > struct rebind
         {
           //! new separable type.
           typedef RegulatedFit< TT_POLICY > type;
         };
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

            //! Computes \a _A and \a _b excluding excluded structures.
            template< class T_MATRIX, class T_VECTOR >
              void other_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const {}
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
            //! alpha for pair regulation
            types::t_real alpha;
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

            //! Adds pair regulation to \a _A.
            template< class T_MATRIX, class T_VECTOR >
              void other_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const;

          protected:
            //! Possible initialization stuff.
            void init( const BaseFit::t_Clusters &_clusters );

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
} // namespace LaDa

#include "fit.impl.h"

#endif 

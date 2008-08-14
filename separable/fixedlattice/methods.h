//
//  Version: $Id$
//
#ifndef _CE_METHODS_H_
#define _CE_METHODS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <crystal/structure.h>

namespace CE
{
  namespace Method
  {
    //! \cond
    namespace Policy
    {
      class Nothing
      {
        public:
          //! Starts policy.
          void start() const {};
          //! Continues policy.
          template< class T_COLLAPSE, class T_MINIMIZER >
            bool go( const T_COLLAPSE &,
                     const T_MINIMIZER &,
                     types::t_int  ) const
              { return false; };
          //! Ends policy.
          template< class T_COLLAPSE, class T_MINIMIZER >
            opt::ErrorTuple end( T_COLLAPSE &_collapse,
                                 T_MINIMIZER &_minimizer, 
                                 types::t_int ) const
            { return _minimizer( _collapse.efficients(), _collapse ); };
      };
      template< class T_COEFFICIENTS >
      class BestOf
      {
        public:
          //! Type of the saved object.
          typedef T_COEFFICIENTS t_Coefficients;
          //! Which best to look at.
          types::t_unsigned which;
          //! How many restarts.
          types::t_unsigned restarts;
          //! How random is random.
          types::t_real howrandom;
          //! Starts policy.
          void start() { nbrestarts = 0; }
          //! Continues policy.
          template< class T_COLLAPSE, class T_MINIMIZER >
            bool go( T_COLLAPSE &,
                     T_MINIMIZER &,
                     types::t_int  );
          //! Ends policy.
          template< class T_COLLAPSE, class T_MINIMIZER >
            opt::ErrorTuple end( T_COLLAPSE &_col,
                                 T_MINIMIZER &, 
                                 types::t_int ) const;

        protected:
          //! Current number of restarts.
          types::t_unsigned nbrestarts;
          //! Current number of restarts.
          opt::ErrorTuple best;
          //! Saved best trial.
          t_Coefficients object;
          //! Save normalization.
          std::vector< typename t_Coefficients::value_type > norms;
      };
    }
    //! \endcond
    //! Class for fitting separables and mixed approaches.
    template< class T_POLICY = Policy::Nothing >
    class Fit
    {
      public:
        //! Type of the policy.
        typedef T_POLICY t_Policy;
        //! Type of the structure set.
        typedef std::vector<Crystal::Structure> t_Structures;
        //! Verbosity.
        types::t_int verbosity;
        //! The policy.
        t_Policy policy;

        //! Constructor.
        Fit( t_Structures &_structures ) : structures_( &_structures), verbosity(false) {}

        template< class T_COLLAPSE, class T_MINIMIZER >
          opt::ErrorTuple operator()( T_COLLAPSE &_collapse, T_MINIMIZER &_minimizer )
          { 
            __TRYBEGIN
            policy.start();
            while( policy.go( _collapse, _minimizer, verbosity - 1 ) );
            opt::ErrorTuple errors( policy.end( _collapse, _minimizer, verbosity - 1 ) );
            if( verbosity >= 1 ) return check_all( _collapse, structures(), verbosity >= 2 );
            return errors;
            __TRYEND(,"Error in CE::Methods::fit().\n" )
          }

        //! Reference to the structure set.
        t_Structures& structures() { return *structures_; }
        //! Constant reference to the structure set.
        const t_Structures& structures() const { return *structures_; }

      protected:
        //! Pointer to the structure set.
        t_Structures *structures_;
    };

    //! \brief Fits a sum of separable functions using specief collapse functor
    //!        and minimizer.
    template< class T_COLLAPSE, class T_MINIMIZER, class T_STRUCTURES >
      opt::ErrorTuple fit( T_COLLAPSE &_collapse,
                           const T_MINIMIZER &_min,
                           const T_STRUCTURES &_strs,
                           bool _verbose = false );
    //! \brief Performs leave-one-out for a sum of separable functions using
    //!        specief collapse functor and minimizer.
    template< class T_COLLAPSE, class T_FIT, class T_MINIMIZER >
      opt::t_ErrorPair leave_one_out( T_COLLAPSE &_collapse,
                                      T_FIT &_fit,
                                      const T_MINIMIZER &_min,
                                      types::t_int _verbosity );
    //! Checks error for one structure.
    template< class T_COLLAPSE >
      opt::ErrorTuple check_one( const T_COLLAPSE &_collapse,
                                 const Crystal::Structure &_structure,
                                 size_t _n, bool _verbose = false );
    //! Checks error for all structure.
    template< class T_COLLAPSE, class T_STRUCTURES >
      opt::ErrorTuple check_all( const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_strs,
                                 bool _verbose = false );
  } // end of namespace Methods.
} // end of namespace CE.

namespace Fitting
{
  //! \brief Makes an alternating linear least-square method from templated 1d
  //!        least-square method.
  //! \tparam T_SOLVER A linear solver type.
  template< class T_SOLVER >
  class AlternatingLeastSquare
  {
    public:
      //! Type of the Matrix for collapsed 1d problem.
      typedef T_SOLVER t_Solver;
      //! Type of the Matrix for collapsed 1d problem.
      typedef boost::numeric::ublas::matrix<types::t_real> t_Matrix;
      //! Type of the Vector for collapsed 1d-problem.
      typedef boost::numeric::ublas::vector<types::t_real> t_Vector;
      //! Type of the input vectors. One for each collapsible dimension.
      typedef std::vector< t_Vector > t_Vectors;

    public:
      //! Maximum number of iterations.
      types::t_int itermax;
      //! Convergence criteria.
      types::t_real tolerance;
      //! The linear solver
      t_Solver linear_solver;
      //! Verbosity
      bool verbose;

    public:
      //! Constructor.
      AlternatingLeastSquare() : itermax(20), tolerance( 1e-4 ), verbose( false )
        { linear_solver.itermax = itermax; linear_solver.tolerance = tolerance;  }
      //! Copy Constructor.
      AlternatingLeastSquare   ( const AlternatingLeastSquare &_c ) 
                             : itermax(_c.itermax), tolerance( _c.tolerance ),
                               linear_solver( _c.linear_solver ),
                               verbose( _c.verbose ) {}
      //! Destructor
      ~AlternatingLeastSquare() {}
      //! \brief Perform Alternating Linear Least-Square Fit.
      template< class T_COLLAPSE >
        opt::ErrorTuple operator()( typename T_COLLAPSE :: t_Matrix& _solution,
                                    T_COLLAPSE &_collapse ) const;
      //! Loads parameters from XML element
      void Load( const TiXmlElement &_node );
  };

}
#include "methods.impl.h"

#endif

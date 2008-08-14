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
    //! \brief Fits a sum of separable functions using specief collapse functor
    //!        and minimizer.
    template< class T_COLLAPSE, class T_MINIMIZER, class T_STRUCTURES >
      opt::ErrorTuple fit( T_COLLAPSE &_collapse,
                           const T_MINIMIZER &_min,
                           const T_STRUCTURES &_strs,
                           bool _verbose = false );
    //! \brief Performs leave-one-out for a sum of separable functions using
    //!        specief collapse functor and minimizer.
    template< class T_COLLAPSE, class T_MINIMIZER, class T_STRUCTURES >
      opt::t_ErrorPair leave_one_out( T_COLLAPSE &_collapse,
                                      const T_MINIMIZER &_min,
                                      const T_STRUCTURES &_strs,
                                      types::t_int _verbosity = 1 );
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

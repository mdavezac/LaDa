//
//  Version: $Id$
//
#ifndef _OPT_ALLSQ_H_
#define _OPT_ALLSQ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iomanip>
#include <vector>


#include <opt/types.h>
#include <opt/debug.h>
#include <tinyxml/tinyxml.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
// #include <gsl/gsl_linalg.h>
// #include <opt/gsl.h>


namespace Fitting
{
  //! \brief Makes an alternating linear least-square method from templated 1d
  //!        least-square method.
  //! \tparam T_SOLVER A linear solver type.
  template< class T_SOLVER >
  class Allsq
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
      Allsq() : itermax(20), tolerance( 1e-4 ), verbose( false )
        { linear_solver.itermax = itermax; linear_solver.tolerance = tolerance;  }
      //! Copy Constructor.
      Allsq   ( const Allsq &_c ) 
            : itermax(_c.itermax), tolerance( _c.tolerance ),
              linear_solver( _c.linear_solver ), verbose( _c.verbose ),
              targets( _c.targets ) {}
      //! Destructor
      ~Allsq() {}
      //! \brief Pre-minimizer stuff.
      //! \param _v[in] is a vector with as many elements as there observed
      //!               data points.
      template< class T_CONTAINER > void init_targets( T_CONTAINER &_v );
      //! \brief Perform Alternating Linear Least-Square Fit.
      //! \params _solution should be a vector of vectors. The outer vectors
      //!         have the dimension of the problem. The inner vectors should
      //!         be of appropriate type for the chosen least-square fit
      //!         method. On output, \a _solution contains the solution.
      //!         On input _solution is used as the starting guess.
      //!         \a _solution is an  std::vector of "something", where
      //!         "something" is a type appropriate for the least-square-fit
      //!         method specified by \a T_LLSQ. 
      //! \tparam T_COLLAPSED is a function type or a functor providing
      //!         void (Collapsed::t_Matrix&, types::t_unsigned, Collapsed::t_Vectors& ).
      //!         - The first argument is a matrix suitable for a 1d
      //!           least-square-fit as required by the templated
      //!           Collapsed::t_Llsq method. 
      //!         - The second argument is dimension which will next undergo
      //!           the least-square fit.
      //!         - The last argument is a reference to the \a _solution, eg
      //!           the current solution vector.
      //!         .
      template< class T_COLLAPSE >
        types::t_real operator()( t_Vectors& _solution, T_COLLAPSE* collapse );
      //! Loads parameters from XML element
      void Load( const TiXmlElement &_node );

    protected:
      //! The target values.
      t_Vector targets;
  };

  template< class T_SOLVER > template< class T_CONTAINER >
    inline void Allsq<T_SOLVER> :: init_targets( T_CONTAINER &_v )
      { 
        targets.resize( _v.size() ); 
        std::copy( _v.begin(), _v.end(), targets.begin() );
      }; 

  template< class T_SOLVER > template< class T_COLLAPSE >
    types::t_real Allsq<T_SOLVER> :: operator()( t_Vectors& _solution,
                                                 T_COLLAPSE* collapse  )
    {

      try
      {
        types::t_real convergence( 1e1 * tolerance );
        if( verbose ) std::cout << "Starting Alternating-least-square fit.\n";
        types::t_unsigned iter = 0;
        types::t_int D( _solution.size() );
        t_Matrix A;
        t_Vector b;
        typename t_Vectors :: iterator i_sol;
        typename t_Vectors :: iterator i_sol_end = _solution.end();
        types::t_real oldvalue;
        do
        {
          types::t_unsigned dim(0);
          i_sol = _solution.begin();
          for(; i_sol != i_sol_end; ++i_sol, ++dim )
          {
            (*collapse)( b, A, dim, targets, _solution );
          // for( size_t i(0); i < A.size1(); ++i )
          // {
          //   for( size_t j(0); j < A.size2(); ++j )
          //     std::cout << std::scientific 
          //               << std::setw(8) 
          //               << std::setprecision(4)
          //               << A(i,j) << " ";
          //   std::cout << "      " 
          //             << std::scientific 
          //             << std::setw(8) 
          //             << std::setprecision(4)
          //             << b(i) << "\n";
          // }
          // std::cout << "\n";
            linear_solver( A, *i_sol, b );
          }
          ++iter;

          if( tolerance > 0e0  )
          {
            types::t_real newvalue = collapse->evaluate( _solution, targets );
            convergence = oldvalue - newvalue;
            oldvalue = newvalue;
            if( verbose )
            {
              if( iter == 1 )
                std::cout << "\n  Allsq iter: " << iter
                          << "  evaluation: " << oldvalue << "\n";
              else
                std::cout << "\n  Allsq iter: " << iter
                          << "  evaluation: " << newvalue 
                          << "  convergence: " << convergence << "\n";
            } // end of if( verbose )
            if( iter > 1 and std::abs(convergence) < tolerance ) return newvalue; 
          } // end of if( tolerance > 0e0 )
        }
        while( iter < itermax or itermax == 0 );

        return collapse->evaluate( _solution, targets );
      }
      __CATCHCODE(, "Error encountered in Alternating-least square fit.\n" )
    }

    template< class T_SOLVER > 
      void Allsq<T_SOLVER> :: Load( const TiXmlElement &_node )
      {
        std::string name = _node.Value();
        const TiXmlElement *parent = &_node;
        if( name.compare( "Allsq" ) ) 
         parent = _node.FirstChildElement( "Allsq" );
        __DOASSERT( not parent, "Could not find Allsq tag in input.\n" )
        if( parent->Attribute( "tolerance" ) )
          parent->Attribute( "tolerance", &tolerance );
        if( parent->Attribute( "itermax" ) )
          parent->Attribute( "itermax", &itermax );
      }

}
#endif

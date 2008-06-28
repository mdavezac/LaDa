//
//  Version: $Id$
//
#ifndef _OPT_ALLSQ_H_
#define _OPT_ALLSQ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <opt/types.h>
#include <opt/debug.h>
#include <tinyxml/tinyxml.h>
#include <gsl/gsl_linalg.h>
#include <opt/gsl.h>

namespace Fitting
{
  //! \brief Makes an alternating linear least-square method from templated 1d
  //!        least-square method.
  //! \tparam T_LSQ The 1d least-square method to alternate for.
  class Allsq
  {
    public:
      //! Type of the Matrix for collapsed 1d problem.
      typedef std::vector<types::t_real> t_Matrix;
      //! Type of the Vector for collapsed 1d-problem.
      typedef std::vector<types::t_real> t_Vector;
      //! Type of the input vectors. One for each collapsible dimension.
      typedef std::vector< t_Vector > t_Vectors;

    public:
      //! Maximum number of iterations.
      types::t_int itermax;
      //! Convergence criteria.
      types::t_real tolerance;

    public:
      //! Constructor.
      Allsq() : itermax(20), tolerance( 1e-4 ) {}
      //! Destructor
      ~Allsq() {}
      //! \brief Pre-minimizer stuff.
      //! \param _v[in] is a vector with as many elements as there observed
      //!               data points.
      void init_targets( t_Vector &_v ) { targets = _v; }; 
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

  template< class T_COLLAPSE >
    types::t_real Allsq :: operator()( t_Vectors& _solution, T_COLLAPSE* collapse  )
    {
      std::cout << "Starting point: " << std::endl;
          typename t_Vectors :: iterator i_sol = _solution.begin();
          typename t_Vectors :: iterator i_sol_end = _solution.end();
          for(i_sol = _solution.begin(); i_sol != i_sol_end; ++i_sol)
          {
            std::cout << "      ";
            std::for_each( i_sol->begin(), i_sol->end(),
                           std::cout << boost::lambda::_1 << "   " );
            std::cout << "\n";
          }
          std::cout << "\n\n";
      try
      {
        types::t_unsigned iter = 0;
        types::t_int D( _solution.size() );
        t_Matrix A;
        t_Vector b;
        types::t_real convergence(0);
        t_Matrix copy;
        do
        {
          types::t_unsigned dim(0);
          i_sol = _solution.begin();
          for(convergence = 0e0; i_sol != i_sol_end; ++i_sol, ++dim )
          {
            (*collapse)( b, A, dim, targets, _solution );
           
            copy = A;
            Gsl::Matrix gslA( b.size(), copy );
            std::cout << ( (gsl_matrix*)gslA )->size1 << " "
                      << ( (gsl_matrix*)gslA )->size2 << " "
                      << ( (gsl_matrix*)gslA )->tda << "\n";
            Gsl::Vector gslb( b );
            Gsl::Vector gslx( *i_sol );
            std::cout << "1A: " ;
            std::for_each( A.begin(), A.end(), std::cout << boost::lambda::_1 << " " );
            std::cout << std::endl;;

            gsl_linalg_HH_solve( (gsl_matrix*) gslA, (gsl_vector*) gslb,
                                 (gsl_vector*) gslx );


            types::t_real result( 0 );
            typename t_Matrix::const_iterator i_A = A.begin();
            typename t_Vector :: const_iterator i_b = b.begin();
            typename t_Vector :: const_iterator i_b_end = b.end();
            for(; i_b != i_b_end; ++i_b )
            {
              types::t_real col(0);
              typename t_Vectors :: value_type
                                 :: const_iterator i_var = i_sol->begin();
              typename t_Vectors :: value_type
                                 :: const_iterator i_var_end = i_sol->end();
              for(; i_var != i_var_end; ++i_var, ++i_A )
              {
                __ASSERT( i_A == A.end(), "" );
                col += (*i_A) * (*i_var );
              }
              std::cout << "|" << col << " - " << *i_b << "| = "
                        << std::abs( col - *i_b ) << std::endl;
              result += col;
            } 
            __ASSERT( i_A != A.end(), "" );
            std::cout << "    dim: " << dim << " conv: " << result << std::endl;
            convergence += result;
          }
          for(i_sol = _solution.begin(); i_sol != i_sol_end; ++i_sol)
          {
            std::cout << "      ";
            std::for_each( i_sol->begin(), i_sol->end(),
                           std::cout << boost::lambda::_1 << "   " );
            std::cout << "\n";
          }
          ++iter;
          convergence /= (types::t_real) D;

          std::cout << "iter: " << iter << " conv: " << convergence << std::endl;
        }
        while(     ( convergence > tolerance or tolerance < 0e0 )
               and ( iter < itermax or itermax == 0 ) );

        std::cout << "final conv: " << convergence << std::endl;
        return convergence;
      }
      __CATCHCODE(, "Error encountered in Alternating-least square fit.\n" )
    }

    inline void Allsq :: Load( const TiXmlElement &_node )
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

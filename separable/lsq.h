//
//  Version: $Id$
//
#ifndef _SEPARABLE_LLSQ_H_
#define _SEPARABLE_LLSQ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <opt/types.h>
#include <opt/debug.h>
#include <tinyxml/tinyxml.h>

namespace Fitting
{
  //! \brief Makes an alternating linear least-square method from templated 1d
  //!        least-square method.
  //! \tparam T_LSQ The 1d least-square method to alternate for.
  template< class T_LLSQ >
  class Allsq
  {
    public:
      //! Type of the 1d least-square method.
      typedef T_LLSQ t_LLSQ;
      //! Type of the Matrix for collapsed 1d problem.
      typedef typename T_LLSQ :: t_Matrix t_Matrix;
      //! Type of the Vector for collapsed 1d-problem.
      typedef typename T_LLSQ :: t_Vector t_Vector;
      //! Type of the input vectors. One for each collapsible dimension.
      typedef typename std::vector< t_Vector > t_Vectors;

    public:
      //! Maximum number of iterations.
      types::t_int itermax;
      //! Convergence criteria.
      types::t_real tolerance;
      //! A least-square-method
      t_LLSQ llsq;

    public:
      //! Constructor.
      Allsq() : itermax(20), tolerance( 1e-4 ) {}
      //! Destructor
      ~Allsq() {}
      //! Pre-minimizer stuff.
      void init_b( t_Vector &_v )
        { llsq.init_b( _v ); }; 
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
  };

  template< class T_LLSQ > template< class T_COLLAPSE >
    types::t_real Allsq<T_LLSQ> :: operator()( t_Vectors& _solution, 
                                               T_COLLAPSE* collapse  )
    {
      try
      {
        types::t_unsigned iter = 0;
        types::t_int D( _solution.size() );
        t_Matrix A;
        types::t_real convergence(0);
        do
        {
          typename t_Vectors :: iterator i_sol = _solution.begin();
          typename t_Vectors :: iterator i_sol_end = _solution.end();
          types::t_unsigned dim(0);
          for(convergence = 0e0; i_sol != i_sol_end; ++i_sol, ++dim )
          {
            (*collapse)( A, dim, _solution );
            llsq.init_A( A );
            convergence += llsq( *i_sol );
          }
          ++iter;
          convergence /= (types::t_real) D;
        }
        while(     ( convergence < tolerance or tolerance < 0e0 )
               and ( iter < itermax or itermax > 0 ) );

        return convergence;
      }
      __CATCHCODE(, "Error encountered in Alternating-least square fit.\n" )
    }

  template< class T_LLSQ >
    void Allsq<T_LLSQ> :: Load( const TiXmlElement &_node )
    {
      std::string name = _node.value();
      TiXmlElement *parent = &_node;
      if( name.compare( "ALLSQ" ) ) 
       parent = _node.FirstChildElement( "ALLSQ" );
      __DOASSERT( not parent, "Could not find ALLSQ tag in input.\n" )
      if( parent->Attribute( "tolerance" ) )
        parent->Attribute( "tolerance", &tolerance );
      if( parent->Attribute( "itermax" ) )
        parent->Attribute( "itermax", &itermax );

      llsq.Load( *parent );
    }

}
#endif

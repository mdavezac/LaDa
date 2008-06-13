//
//  Version: $Id$
//
#ifndef _SEPARABLE_LLSQ_H_
#define _SEPARABLE_LLSQ_H_

#include <vector>

#include <opt/types.h>
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
      typedef T_LLSQ :: t_Matrix t_1dMatrix;
      //! Type of the Vector for collapsed 1d-problem.
      typedef T_LLSQ :: t_Vector t_1dVector;
      //! Conmtainer of Matrices for each alternating direction.
      typedef std::vector< t_1dMatrix > t_Matrices;
      //! Conmtainer of Matrices for each alternating direction.
      typedef std::vector< t_1dVector > t_Vectors;

    public:
      //! Maximum number of iterations.
      types::t_int itermax;
      //! Convergence criteria.
      types::t_real tolerance;
      //! A least-square-method
      t_LLSQ llsq;

    public:
      //! Constructor.
      Allsq() : itermax(20), tolerance( 1e-4 )  {}
      //! Destructor
      ~Allsq() {}
      //! Pre-minimizer stuff.
      void init( t_Matrices &_m, t_Vectors &_v ) { init_A(_m); init_b(_v); }; 
      //! Pre-minimizer stuff.
      void init_b( t_1dVector &_v )
        { lsq.init_b( _v ) }; 
      //! Pre-minimizer stuff.
      void init_A( t_Matrices &_m ) { matrices = _m; }; 
      //! \brief Perform Alternating Linear Least-Square Fit.
      //! \params _solution should be a vector of vectors. The outer vectors
      //!         have the dimension of the problem. The inner vectors should
      //!         be of appropriate type for the chosen least-square fit
      //!         method. On output, \a _solution contains the solution.
      //!         On input _solution is used as the starting guess.
      types::t_real operator()( t_Vectors& _solution );
      //! Loads parameters from XML element
      void Load( const TiXmlElement &_node );
      

    protected:
      //! Matrices in each alternating direction.
      t_Matrices matrices;
  };

  template< class T_LLSQ >
    types::t_real Allsq<T_LLSQ> :: operator()( t_Vectors& _solution )
    {
      __DOASSERT(     _solution.size() == matrices.size()
                  and matrices.size() == vectors.size(),
                  "Incoherent sizes of matrices/vectors Ax=b." )
                 
      types::t_unsigned iter = 0;
      types::t_int D( _solution.size() );
      llsq.init( vector );
      do
      {
        types::t_real convergence;
        typename t_Vectors :: iterator i_sol = _solution.begin();
        typename t_Vectors :: iterator i_sol = _solution.end();
        typename t_Matrices :: const_iterator i_A = matrices.end();
        for( convergence = 0e0; i_sol != i_sol_end; ++i_sol, ++i_A )
        {
          llsq.init( *i_A );
          convergence += llsq( *i_sol );
        }
        ++iter;
        convergence /= (types::t_real) D;
      }
      while(     ( convergence < tolerance or tolerance < 0e0 )
             and ( iter < itermax or itermax > 0 ) );
      return convergence;
    }

  template< class T_LLSQ >
    void Allsq<T_LLSQ> :: Load( const TiXmlElement &_node )
    {
      std::string name = _node.value();
      TiXmlElement *parent = &_node;
      if( name.compare( "ALLSQ" ) ) 
       parent = _node.FirstChildElement( "ALLSQ" );
      __DOASSERT( not parent, "Could not find ALLSQ tag in input.\n" )
      if( parent->Attribute( "convergence" ) )
        parent->Attribute( "convergence", &convergence );
      if( parent->Attribute( "itermax" ) )
        parent->Attribute( "itermax", &itermax );

      lsq.Load( *parent );
    }

}
#endif

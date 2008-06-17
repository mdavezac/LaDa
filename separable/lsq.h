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
      //! \brief Whether to update after all dims have been optimized, or in
      //!        between dimensions.
      bool update_between_dims;

    public:
      //! Constructor.
      Allsq() : itermax(20), tolerance( 1e-4 ), update_between_dims(false)  {}
      //! Destructor
      ~Allsq() {}
      //! Pre-minimizer stuff.
      void init_b( t_Vector &_v )
        { lsq.init_b( _v ) }; 
      //! \brief Perform Alternating Linear Least-Square Fit.
      //! \params _solution should be a vector of vectors. The outer vectors
      //!         have the dimension of the problem. The inner vectors should
      //!         be of appropriate type for the chosen least-square fit
      //!         method. On output, \a _solution contains the solution.
      //!         On input _solution is used as the starting guess.
      template< class T_COLLAPSED >
        types::t_real operator()( t_Vectors& _solution, T_COLLAPSE collapse );
      //! Loads parameters from XML element
      void Load( const TiXmlElement &_node );
  };

  template< class T_LLSQ > template< class T_COLLAPSE >
    types::t_real Allsq<T_LLSQ> :: operator()( t_Vectors& _solution, 
                                               T_COLLAPSE collapse  )
    {
      t_Vectors *save(NULL);
      __DOASSERT(     _solution.size() == matrices.size()
                  and matrices.size() == vectors.size(),
                  "Incoherent sizes of matrices/vectors Ax=b." )

      try
      {
        types::t_unsigned iter = 0;
        types::t_int D( _solution.size() );
        t_Matrix A;
        if( not update_between_dims ) save = new t_Vectors( _solution );
        do
        {
          types::t_real convergence;
          typename t_Vectors :: iterator i_sol = _solution.begin();
          typename t_Vectors :: iterator i_sol_end = _solution.end();
          types::t_unsigned dim(0);
          for(convergence = 0e0; i_sol != i_sol_end; ++i_sol, ++dim )
          {
            collapse( A, dim, save ? *save: _solution )
            llsq.init( A );
            convergence += llsq( *i_sol );
          }
          ++iter;
          convergence /= (types::t_real) D;
          if( save ) std::copy( _solution.begin(), _solution.end(), save->begin() );
        }
        while(     ( convergence < tolerance or tolerance < 0e0 )
               and ( iter < itermax or itermax > 0 ) );

        if( save ) delete save;

        return convergence;
      }
      __CATCHCODE( if( save ) delete save; save = NULL;,
                   "Error encountered in Alternating-least square fit." )
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

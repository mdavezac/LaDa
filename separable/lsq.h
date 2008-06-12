//
//  Version: $Id$
//
#ifndef _SEPARABLE_LSQ_H_
#define _SEPARABLE_LSQ_H_

//! Namespace for alternating least-square method.
namespace ALSQ
{
  //! \brief Makes an alternating least-square method from templated 1d
  //!        least-square method.
  //! \tparam T_LSQ The 1d least-square method to alternate for.
  template< class T_LSQ >
  class Base
  {
    public:
      //! Type of the 1d least-square method.
      typedef T_LSQ t_LSQ;
      //! Type of the Matrix for collapsed 1d problem.
      typedef T_LSQ :: t_Matrix t_1dMatrix;
      //! Type of the Vector for collapsed 1d-problem.
      typedef T_LSQ :: t_Vector t_1dVector;
      //! Conmtainer of Matrices for each alternating direction.
      typedef std::vector< t_1dMatrix > t_Matrices;
      //! Conmtainer of Matrices for each alternating direction.
      typedef std::vector< t_1dVector > t_Vectors;

    public:
      //! Maximum number of iterations.
      types::t_int itermax;
      //! Convergence criteria.
      types::t_real epsilon;

    public:
      //! Constructor.
      Base() : itermax(20), epsilon( 1e-4 )  {}
      //! Destructor
      ~Base() {}
      //! Pre-minimizer stuff.
      void init( t_Matrices &_m, t_Vectors &_v ) { matrices = _m; vectors = _v; }; 
      

    protected:
      //! Matrices in each alternating direction.
      t_Matrices matrices;
      //! Vectors in each alternating direction.
      t_Matrices vectors;
      //! A least-square-method
      t_LSQ lsq;
  }
}
#endif

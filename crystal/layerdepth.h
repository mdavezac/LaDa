//
//  Version: $Id$
//

#ifndef _LADA_CRYSATL_LAYERDEPTH_H_
#define _LADA_CRYSATL_LAYERDEPTH_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/debug.h>
#include <math/fuzzy.h>
#include <boost/tuple/tuple.hpp>

#include <math/eigen.h>

namespace LaDa
{

  namespace Crystal 
  {
    //! \brief Strict Weak Ordering functor according to depth along eptiaxial
    //!        direction
    //! \details Two vectors are compared by using the value of their scalar
    //!          product with Depth::matrix.col(0) (mod |matrix.col(0)|). If
    //!          these scalar product are equal (as defined by math::eq()),
    //!          then their scalar product with Depth::a1 are compared. If
    //!          again these are equal, then the scalar porducts with Depth::a2
    //!          are compared and the result return.
    class LayerDepth
    {
      protected:
        math::rMatrix3d matrix_; //!< Ordering directions.
      public:
        //! Copy Constructor
        LayerDepth( const LayerDepth &_c ) : matrix_(_c.matrix_) {}
        //! Constructor and Initializer
        LayerDepth( const math::rVector3d &_vec ) 
        {
          matrix_.col(0) = _vec;
          matrix_.col(1) = math::eq( _vec(0), 0e0 ) ? 
                               ( math::eq( _vec(1), 0e0 ) ? 
                                  math::rVector3d(1, 0, 0):
                                  math::rVector3d(0, _vec(2), -_vec(1) )  ): 
                               math::rVector3d( -_vec(2) -_vec(1), _vec(0), _vec(0) );
          matrix_.col(1).normalize();
          matrix_.col(2) = matrix_.col(0).cross(matrix_.col(1));
          matrix_.col(2).normalize();
        }
        //! Destructor.
        virtual ~LayerDepth() {}
    
        //! Strict weak ordering operator.
        bool operator()( const math::rVector3d& _first,
                         const math::rVector3d& _second ) const;
        //! Returns the depth.
        types::t_real operator()( const math::rVector3d& _first ) const;
    };

  } // namespace Crystal

} // namespace LaDa

#endif

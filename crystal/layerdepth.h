//
//  Version: $Id$
//

#ifndef _LADA_CRYSATL_LAYERDEPTH_H_
#define _LADA_CRYSATL_LAYERDEPTH_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <atat/vectmac.h>

#include <opt/debug.h>
#include <opt/fuzzy.h>
#include <boost/tuple/tuple.hpp>

namespace LaDa
{

  namespace Crystal 
  {
    //! \brief Strict Weak Ordering functor according to depth along eptiaxial
    //!        direction
    //! \details Two vectors are compared by using the value of their scalar
    //!          product with Depth::a0 (mod |a1|). If these scalar product are equal (as
    //!          defined by Fuzzy::eq()), then their
    //!          scalar product with Depth::a1 are compared. If again these are
    //!          equal, then the scalar porducts with Depth::a2 are compared and
    //!          the result return. Depth::a0 is the first column of the matrix
    //!          given in the constructor argument. Depth::a2 is the second
    //!          column, and Depth::a3 the third. 
    class LayerDepth
    {
      protected:
        atat::rVector3d a0; //!< First ordering direction
        atat::rVector3d a1; //!< Second ordering direction
        atat::rVector3d a2; //!< Third ordering direction
        __DODEBUGCODE( bool isset; ) 
    
      public:
        //! Constructor.
        LayerDepth() __DODEBUGCODE( : isset(false) ) {}
        //! Copy Constructor
        LayerDepth   ( const LayerDepth &_c )
                   : a0( _c.a0 ), a1( _c.a1 ), a2( _c.a2 )
                     __DODEBUGCODE( __COMMA__ isset( _c.isset ) ) {}
        //! \brief Constructor and Initializer
        //! \param _mat Depth::a0 is set to this vector.
        //!             Depth::a1 and Depth::a2 are constructed.
        LayerDepth   ( const atat::rVector3d &_vec ) { set( _vec); }
        //! \brief Constructor and Initializer
        //! \param _mat Depth::a0 is set to the (normalized) first column of this matrix,
        //!             Depth::a1 to the second, and Depth::a2 to the third.
        LayerDepth   ( const atat::rMatrix3d &_mat ) { set( _mat); }
        //! \brief Constructor and Initializer
        //! \param _mat Depth::a0 is set to this vector.
        //!             Depth::a1 and Depth::a2 are constructed.
        LayerDepth   ( const boost::tuple<types::t_real, types::t_real, types::t_real> &_vec ) 
         { set( _vec ); }
        //! Destructor.
        virtual ~LayerDepth() {}
        //! Sets reference vectors.
        void set( const atat::rMatrix3d &_mat );
        //! Sets and constructs reference vectors.
        void set( const atat::rVector3d &_mat );
        //! Sets and constructs reference vectors.
        void set( const boost::tuple<types::t_real, types::t_real, types::t_real> &_vec ) 
          { set( atat::rVector3d( _vec.get<0>(), _vec.get<1>(), _vec.get<2>() ) ); }
    
        //! Strict weak ordering operator.
        bool operator()( const atat::rVector3d& _first,
                         const atat::rVector3d& _second ) const;
        //! Returns the depth.
        types::t_real operator()( const atat::rVector3d& _first ) const;
    };

  } // namespace Crystal

} // namespace LaDa

#endif

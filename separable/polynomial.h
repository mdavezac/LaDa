//
//  Version: $Id$
//
#ifndef _SEPARABLE_POLYNOMIAL_H_
#define _SEPARABLE_POLYNOMIAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Create polynomials from degree -21 to degree 21
//! \cond
#     define __POLYFUNCTIONS__
#     include "polynomial.impl.h"
//! \endcond

//! Contains all things separable functions.
namespace Separable
{
  //! A polynomial function with integer power.
  struct Polynomial
  {
    //! Type of the argument.
    typedef types::t_double t_Arg;
    //! Type of the return.
    typedef types::t_double t_Return;
    //! Degree of the polynomial.
    types::t_int degree;
    //! \f$\_x^{\text{degree}}$\f.
    t_Return operator()( t_Arg _x ) const
      { return details::Pow<degree>::pow<t_Arg>( _x ); }
    //! \f$\text{degree}\cdot\_x^{\text{degree}-1}$\f.
    t_Return gradient( t_Arg _x ) const
      { return details::Pow<degree>::gradient<t_Arg>( _x ); }
  };

  //! A basis of polynomials.
  class Polynomials 
  {
    protected:
      //! Type of the container of functions.
      typedef std::vector< Polynomial > t_Container;

    public:
      //! Type of a polynomial function.
      typedef t_Container :: value_type value_type;
      //! Type of an iterator for the basis.
      typedef t_Container :: iterator iterator;
      //! Type of a constant iterator for the basis.
      typedef t_Container :: const_iterator const_iterator;

      //! Constructor. Specifies extent of basis.
      Polynomial   ( types::t_int _min, types::t_int _max )
                 : min( _min ), max( _max ), container( _max - _min )
        { set_range( _min, _max ); }
      //! Destructor.
      ~Polynomial() {};
      //! Returns range of basis
      std::pair< types::t_int, types::t_int > get_range() const 
        { return std::pair< types::t_int, types::t_int >( min, max ); }
      //! Sets range of basis
      void set_range( std::pair< types::t_int, types::t_int > _p )
        { set_range( _p.first, _p.second ); }
      //! Sets range of basis
      void set_range( types::t_int _min, types::t_int _max );

    protected:
      //! Container of polynomials.
      t_Container container;
      //! Minimum degree in basis.
      types::t_int min;
      //! Maximum degree in basis.
      types::t_int max;
  };

} // end of Separable namespace

#endif

//
//  Version: $Id$
//
#ifndef _HARMONICS_H_
#define _HARMONICS_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <math.h>

#include <tinyxml/tinyxml.h>


#include <opt/types.h>
#include <opt/traits.h>
#include <atat/vectmac.h>

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace Ising_CE 
{

  //! \brief Creates a %function through the linear interpolation of a set of points.
  //! \details Points are defined in 2D real space. This class evaluates the
  //!          value of the function by finding the two closest points on \e x
  //!          axis and performing a linear interpolation. It can also compute
  //!          the gradient in a similar way.
  class Linear_Interpolator 
  {
    //! \brief Defines a 2-dimensional point.
    //! \details Defines an x and y member, as well as comparison functions for
    //!          easily retriving the closes instances to a point.
    struct Point
    {
      //! Coordinate along the \e x axis.
      types::t_real x;
      //! Coordinate along the \e y axis.
      types::t_real y;

      //! Constructor
      Point() {}
      //! Copy Constructor
      Point(const Point &_p) : x(_p.x), y(_p.y)  {}
      //! Constructor and initializer
      Point(const types::t_real _x, const types::t_real _y) : x(_x), y(_y)  {}
      //! Compares (Fuzzy math) the \e x axis coordinate with \a _x
      bool x_greater( const types::t_real _x ) const
        { return opt::Fuzzy<types::t_real>::greater(x, _x); }
      //! Compares (Fuzzy math) the \e x axis coordinate with \a _x
      bool x_lesser( const types::t_real _x ) const 
        { return opt::Fuzzy<types::t_real>::less(x, _x); }
    };

#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<Linear_Interpolator> ( Linear_Interpolator& );
    //! \endcond
#endif
    

    protected:
      //! A collection of Points between which to interpolate.
      std::vector<Point> points;

    public:
      //! Constructor.
      Linear_Interpolator(){};
      //! Copy Constructor.
      Linear_Interpolator(const Linear_Interpolator &_c) : points(_c.points) {}

      //! Adds a point to the interpolator.
      void add_point( const types::t_real _x, const types::t_real _y);

      //! Returns the interpolated value at \a _x.
      types::t_real evaluate( const types::t_real _x ) const;
      //! Returns the interpolated gradient at \a _x.
      types::t_real evaluate_gradient( const types::t_real _x ) const;
      //! Returns the interpolated value and the computes the gradient at \a _x.
      types::t_real evaluate_with_gradient( const types::t_real _x, types::t_real &_grad ) const;

      //! Loads the set of points to interpolate from XML
      bool Load (const TiXmlElement &_element);

      //! Cleats the set of interpolation point.
      void clear() { points.clear(); };
  };


  /** \brief Defines a cubic harmonic %function of rank 0, 1, 2, or 3.
   *  \details This class can also compute something  closer to what the
   *           constituent strain requires: the product between an interpolated
   *           %function of the concentration, a gaussian attenuation of the
   *           norm of the reciprocal-space vector, and the cubic haromic proper.
   *           With \f$\hat{k}\f$ the normalized reciprocal-space vector, the
   *           expression for the cubic-harmonic proper are:
   *    - Rank 0: \f$ \sqrt{\frac{1}{4\pi}} \f$
   *    - Rank 1: \f$ -\sqrt{\frac{21}{\pi}}\left[ 1 - 5
   *                  \left(\hat{k}_x^2\hat{k}_y^2 + \hat{k}_x^2 \hat{k}_z^2 +
   *                  \hat{k}_y^2 \hat{k}_z^2\right)\right] \f$
   *    - Rank 2:  \f$ \frac{1}{4}\sqrt{\frac{13}{8\pi}}
   *                   \left[  7 \left(   \hat{k}_x^6 + \hat{k}_y^6 + \hat{k}_z^6 
   *                                     + 30 \hat{k}_x^2\hat{k}_y^2\hat{k}_z^2 
   *                             \right) - 5
   *                   \right] \f$
   *    - Rank 3: \f$ \frac{1}{8} \sqrt{\frac{561}{4\pi}}
   *                  \left[   \hat{k}_x^8 + \hat{k}_y^8 + \hat{k}_z^8 
   *                         - 14 \left(   \hat{k}_x^6\hat{k}_x^2 
   *                                     + \hat{k}_x^6\hat{k}_z^2 
   *                                     + \hat{k}_y^6\hat{k}_z^2 
   *                                     + \hat{k}_x^2\hat{k}_y^6 
   *                                     + \hat{k}_x^2\hat{k}_z^6 
   *                                     + \hat{k}_y^2\hat{k}_z^6 
   *                              \right)
   *                         + 35 \left(   \hat{k}_x^4 \hat{k}_y^4
   *                                     + \hat{k}_x^4 \hat{k}_z^4
   *                                     + \hat{k}_y^4 \hat{k}_z^4
   *                              \right) 
   *                  \right] \f$
   *    .
   *           For more information on harmonics, constituent strain, or the
   *           %Cluster Formalism, you can start here:  <A
   *           HREF="http://dx.doi.org/10.1103/PhysRevB.46.12587"> David B.
   *           Laks, \e et \e al. PRB \b 46, 12587-12605 (1992) </A>.
   */                                    
  class Harmonic 
  {
#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<Harmonic> ( Harmonic& );
    //! \endcond
#endif

    protected:
      //! Interpolation %function for the harmonic
      Linear_Interpolator interpolation;
      //! Rank of the harmonic
      types::t_unsigned rank;
      //! Attenuation in reciprocal space of the harmonic
      static types::t_real attenuation;

    public:
      //! Constructor
      Harmonic() {};
      //! Copy Constructor
      Harmonic   ( const Harmonic &_h )
               : interpolation( _h.interpolation ),
                 rank( _h.rank ) {}

    public:

      //! \brief Computes the product of the interpolation at \a _x, of a
      //!    traitsgaussian atternuation at \a _k , and of the cubic harmonic at \a _k.
      //! \param [in] _x concentration.
      //! \param [in] _k reciprocal-space vector.
      types::t_real evaluate(const types::t_real _x, const atat::rVector3d &_k) const;
      //! \brief Returns the interpolated value at \a _x.
      //! \param [in] _x concentration.
      types::t_real evaluate(const types::t_real _x) const
        { return (  interpolation.evaluate(_x) ); }
      //! \brief Computes the product of a gaussian attenuation at \a _k and of
      //!        the cubic harmonic at \a _k.
      //! \param [in] _k reciprocal-space vector.
      types::t_real evaluate(const atat::rVector3d &_k) const
        { return (  exp( -norm2(_k) * attenuation ) * evaluate_harmonic( _k ) ); }
      //! \brief returns the gradient of the product of the interpolation at \a _x, of a
      //!        gaussian atternuation at \a _k , and of the cubic harmonic at \a _k.
      //! \param [in] _x concentration.
      //! \param [in] _k reciprocal-space vector.
      types::t_real evaluate_gradient(const types::t_real _x, const atat::rVector3d &_k) const;
      //! \brief Returns the interpolated gradient at \a _x.
      //! \param [in] _x concentration.
      types::t_real evaluate_gradient(const types::t_real _x) const
        { return ( interpolation.evaluate_gradient(_x) ); }
      //! \brief returns the value and computes the gradient of the product of
      //!        the interpolation at \a _x, of a gaussian atternuation at \a
      //!        _k , and of the cubic harmonic at \a _k.
      //! \param [in] _x concentration.
      //! \param [in] _k reciprocal-space vector.
      //! \param [in, out] _grad stores the computed gradient.
      types::t_real evaluate_with_gradient(const types::t_real _x,
                                           const atat::rVector3d &_k, 
                                           types::t_real &_grad) const;
      //! \brief Returns the interpolated value and computes the gradient at \a _x.
      //! \param [in] _x concentration.
      //! \param [in, out] _grad stores the computed gradient.
      types::t_real evaluate_with_gradient(const types::t_real _x,
                                           types::t_real &_grad) const
        { return interpolation.evaluate_with_gradient(_x, _grad); }

      //! \brief Returns the value \a _k of the cubic harmonic.
      //! \param [in] _k reciprocal-space vector.
      types::t_real evaluate_harmonic(const atat::rVector3d &_k) const;
      
      //! Loads the cubic harmonic and the interpolation from XML.
      bool Load(const TiXmlElement &_element);

      //! Sets the attenuation.
      static void set_attenuation( const types::t_real _a ) 
        { attenuation = (_a == 0) ? 0 : 1.0 / (_a*_a); }

      //! Clears the interpolation.
      void clear() { interpolation.clear(); };
  };

  inline types::t_real Harmonic :: evaluate(const types::t_real _x,
                                            const atat::rVector3d &_k) const
  {
    return (   interpolation.evaluate(_x) 
             * exp( -norm2(_k) * attenuation )
             * evaluate_harmonic( _k ) );
  }
  inline types::t_real Harmonic :: evaluate_gradient(const types::t_real _x,
                                                     const atat::rVector3d &_k) const
  {
    return (   interpolation.evaluate_gradient(_x) 
             * exp( -norm2(_k) *  attenuation )
             * evaluate_harmonic( _k ) );
  }
  inline types::t_real Harmonic :: evaluate_with_gradient(const types::t_real _x, 
                                                   const atat::rVector3d &_k,
                                                   types::t_real &_grad) const
  {
    types::t_real factor = exp( -norm2(_k) * attenuation )
                    * evaluate_harmonic( _k );
    types::t_real result =   interpolation.evaluate_with_gradient(_x, _grad) 
                    * factor ;
    _grad *= factor;
    return result;
  }

} // namespace Ising_CE 

#ifdef _MPI

namespace mpi
{
  /** \ingroup MPI
   *  \brief serializes an Ising_CE::Linear_Interpolator
   *  \details This includes serializing all points ant their x and y members. */
  template<>
  bool BroadCast :: serialize<Ising_CE::Linear_Interpolator>
                             ( Ising_CE::Linear_Interpolator &_l);
  /** \ingroup MPI
   *  \brief serializes an Ising_CE::Harmonic
   *  \details This includes serializing the rank, the attenuation, and the
   *           interpolation. */
  template<>
  bool BroadCast :: serialize<Ising_CE::Harmonic>( Ising_CE::Harmonic &_h );
}

#endif
#endif // _HARMONICS_H_

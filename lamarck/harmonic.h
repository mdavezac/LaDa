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
#include <opt/fuzzy.h>
#include <atat/vectmac.h>
#include <mpi/mpi_object.h>

#include <mpi/mpi_object.h>

namespace Ising_CE 
{
  namespace ConstituentStrain
  {
    //! \brief Constains all things Harmonics.
    //! \details Defines a linear interpolator which interpolates between
    //!          concentration values. Also defines the harmonics via a
    //!          curiously recurring template mechanism. The base class
    //!          Harmonic::Base defines most all functions necessary for the
    //!          constituent strain minus the evaluation of the harmonic
    //!          proper. That is left to the derived classes. 
    //!          For more information on harmonics, constituent strain, or the
    //!          %Cluster Expansion Formalism, you can start here:  <A
    //!          HREF="http://dx.doi.org/10.1103/PhysRevB.46.12587"> David B.
    //!          Laks, \e et \e al. PRB \b 46, 12587-12605 (1992) </A>.
    namespace Harmonic
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
            { return Fuzzy::gt(x, _x); }
          //! Compares (Fuzzy math) the \e x axis coordinate with \a _x
          bool x_lesser( const types::t_real _x ) const 
            { return Fuzzy::le(x, _x); }
          //! Serializes a point.
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
            { _ar & x; _ar & y; } 
        };
    
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
          types::t_real evaluate_with_gradient( const types::t_real _x,
                                                types::t_real &_grad ) const;
    
          //! Loads the set of points to interpolate from XML
          bool Load (const TiXmlElement &_element);
    
          //! Cleats the set of interpolation point.
          void clear() { points.clear(); };
          //! Serializes a linear interpolation.
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
            { _ar & points; } 
      };
    
      //! \brief Base class for harmonics.
      //! \details Allows static "virtual" functions \e via curriously recurring templates.
      //!          This class can also compute something  closer to what the
      //!          constituent strain requires: the product between an interpolated
      //!          %function of the concentration, a gaussian attenuation of the
      //!          norm of the reciprocal-space vector, but not the cubic haromic proper.
      //!          That particular implementation is left to the derived classes.
      //!          Derived classes must contain the following members
      //!             - types::t_real operator()(const atat::rVector3d& _k ) const
      //!               returning the harmonic in direction \e _k.
      //!             - const static std::string type describing with one word
      //!               the type of harmonic. This variable is used when
      //!               loading harmonics from XML.
      //!             - const static types::t_int maxrank contains the maximum
      //!               rank implemented for that harmonic. This rank is
      //!               checked for when loading from XML.
      //!             .
      template< class T_DERIVED >
      class Base   
      {
        public:
          //! type of the derived class
          typedef T_DERIVED t_Derived;
        protected:
          //! Interpolation %function for the harmonic
          Linear_Interpolator interpolation;
          //! Rank of the harmonic
          types::t_unsigned rank;
          //! Attenuation in reciprocal space of the harmonic
          static types::t_real attenuation;
    
        public:
          //! Constructor
          Base() {};
          //! Copy Constructor
          Base   ( const Base<t_Derived> &_h )
               : interpolation( _h.interpolation ), rank( _h.rank ) {}
    
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
          types::t_real evaluate(const atat::rVector3d &_k) const;
          //! \brief returns the gradient of the product of the interpolation at \a _x, of a
          //!        gaussian atternuation at \a _k , and of the cubic harmonic at \a _k.
          //! \param [in] _x concentration.
          //! \param [in] _k reciprocal-space vector.
          types::t_real evaluate_gradient(const types::t_real _x,
                                          const atat::rVector3d &_k) const;
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
    
          //! Loads the cubic harmonic and the interpolation from XML.
          bool Load(const TiXmlElement &_element);
    
          //! Sets the attenuation.
          static void set_attenuation( const types::t_real _a ) 
            { attenuation = (_a == 0) ? 0 : 1.0 / (_a*_a); }
    
          //! Clears the interpolation.
          void clear() { interpolation.clear(); };
          //! Serializes a harmonic.
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
            { _ar & interpolation; _ar & rank; _ar & attenuation; } 
        protected:
          //! Ugly hack to deal with inconsistent tetragonal lattice used in CE@nrel.
          void transform_k( const atat::rVector3d &_in, atat::rVector3d &_out ) const 
            { _out = _in; }

      };
      template< class T_DERIVED > types::t_real Base<T_DERIVED> :: attenuation = 1000000.0;
    
      /** \brief Defines a cubic harmonic %function of rank 0, 1, 2, or 3.
       *  \details With \f$\hat{k}\f$ the normalized reciprocal-space vector, the
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
       */                                    
      class Cubic : public Base< Cubic >
      {
        friend class Base< Cubic >;
        public:
          //! Constructor
          Cubic() {};
          //! Copy Constructor
          Cubic   ( const Cubic &_h )
                : Base<Cubic>( _h ) {}
    
        public:
          //! \brief Returns the value \a _k of the cubic harmonic.
          //! \param [in] _k reciprocal-space vector.
          types::t_real operator()(const atat::rVector3d &_k) const;
          
        public:
          //! Names the Type of this harmonic
          const static std::string type;
          //! Maximum rank of the cubic harmonic
          const static types::t_int maxrank = 3;
      };

      //! Defines a tetragonal harmonic %function of rank 0-9.
      class Tetragonal : public Base< Tetragonal >
      {
        friend class Base< Tetragonal >;
        public:
          //! Constructor
          Tetragonal() {};
          //! Copy Constructor
          Tetragonal   ( const Tetragonal &_h ) 
                     : Base<Tetragonal>( _h ) {}
    
        public:
          //! \brief Returns the value \a _k of the cubic harmonic.
          //! \param [in] _k reciprocal-space vector.
          types::t_real operator()(const atat::rVector3d &_k) const;
          
        public:
          //! Names the Type of this harmonic
          const static std::string type;
          //! Maximum rank of the cubic harmonic
          const static types::t_int maxrank = 9;
        protected:
          //! Ugly hack to deal with inconsistent tetragonal lattice used in CE@nrel.
          void transform_k( const atat::rVector3d &_in, atat::rVector3d &_out ) const 
            { _out[0] = _in[0], _out[1] = _in[1]; _out[2] = _in[2] * 1.2;  }
      };
    
    } // namespace Harmonic
  
  } // namespace ConstituentStrain


} // namespace Ising_CE 

#include "harmonic.impl.h"
#endif // _HARMONICS_H_

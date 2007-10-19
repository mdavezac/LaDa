//
//  Version: $Id$
//
#ifndef _CONCENETRATION_H_
#define _CONCENETRATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#include "lamarck/structure.h"
#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif


/** \ingroup Genetic
 * @{ */
//! \brief  Creates a function \f$x(y)\f$ for load matching in ternaries and !
//! quaternaries (see <A HREF="http://dx.doi.org/10.1063/1.2010621"> Rita
//! Magri, Alex Zunger, H Kroemer, J. App. Phys. <STRONG>98</STRONG>, 43701 (2005)</A>)
//! \details For a quaternary
//! \f$\mathrm{A}_x\mathrm{B}_y\mathrm{C}_{1-x}\mathrm{D}_{1-y}\f$ this equation
//! gives us \f$x(y)\f$ such that the inplane stress is (near) zero. The atoms A
//! and C are the first and second atoms of the first site of a lattice, whereas
//! B and D are the first and second atoms of the second site. In other words, we
//! expect the following input
//!   \code
// <Lattice>
//   <column x="0.0" y="0.5" z="0.5" />
//   <column x="0.5" y="0.0" z="0.5" />
//   <column x="0.5" y="0.5" z="0.0" />
//   <site x="0.25" y="0.25" z="0.25">
//       <atom type="A" />
//       <atom type="C" />
//   </site>
//   <site x="0.0" y="0.0" z="0.0">
//       <atom type="B" />
//       <atom type="D" />
//   </site>
// </Lattice>
// <Functional type="Concentration" a="0.001" b="0.648" c="0.239" />
//!   \endcode
//! The function is a second order polynomial given  by \f$x(y) = a\cdot y^2 +
//! b\cdot y + c\f$, where the parmaters are given on input as in the XML example
//! above.
  class X_vs_y
  {
#ifdef _MPI
    //! allows the serialization of an X_vs_y object
    friend bool mpi::BroadCast::serialize<X_vs_y>( X_vs_y & );
#endif
    protected:
      types::t_real a; //!< Second order parameter, as in \f$a\cdot y^2\f$
      types::t_real b; //!< First order parameter, as in \f$b\cdot y\f$
      types::t_real c; //!< zero order parameter
      types::t_real x0; //!< Set x concentration
      types::t_real y0; //!< Set y concentration
    public:
      bool single_c; //!< true if concentration is fixed

    public:
      //! Constructor
      X_vs_y() : a(0), b(0), c(0), x0(0), y0(0), single_c(false) {}
      //! Copy Constructor
      X_vs_y( const X_vs_y &_c) : a(_c.a), b(_c.b), c(_c.c), 
                                  x0(_c.x0), y0(_c.y0), single_c(_c.single_c) {}
      //! Destructor
      ~X_vs_y () {}

      //! \brief Loads parameters from XML
      //! \details The input is expected to appear in the following form if the
      //! concentration is not set.
      //! \code
      // <Functional type="Concentration" a="0.001" b="0.648" c="0.239" />
      //! \endcode 
      //! Additionally, to fix the concentration, an x0, or a y0 attribute can
      //! be added. If both an x0 and a y0 attribute are found, then a, b, c
      //! are irrelevant.
      bool Load( const TiXmlElement &_node );
      //! \brief Returns x with respect to input _y
      //! \details returns X_vs_y::x0 if the concentration is fixed.
      types::t_real get_x( types::t_real _y );
      //! \brief sets a fixed concentration
//     void set_xy( types::t_real _x, types::t_real _y )
//       { x0 = _x; y0 = _y; single_c = true; }
      //! returns X_vs_y::y0
      types::t_real get_y() { return y0; }
      //! returns X_vs_y::x0
      types::t_real get_x() { return x0; }
      //! \brief Returns y with respect to input _x
      //! \details returns X_vs_y::y0 if the concentration is fixed.
      types::t_real get_y( types::t_real _x );
      //! Returns true if some y can be computed with respect to input _x
      bool can_inverse( types::t_real _x );
      //! Returns true if the concentration is fixed
      bool is_single_c () const { return single_c; }

      //! Returns a string describing this class
      std::string print() const;
  };

  inline types::t_real X_vs_y :: get_x( types::t_real _y ) 
  {
    if ( single_c ) return x0;
    return c + b * _y + a * _y * _y; 
  }
  
  inline types::t_real X_vs_y :: get_y( types::t_real _x )
  {
    if ( single_c ) return y0;
    if ( std::abs ( a ) < types::tolerance )
      return ( _x - c ) / b;
   
    types::t_real det = b*b - 4.0 * (c-_x) * a; 
    if ( det < 0 )
    {
      std::cerr << "Error when using Concentration::get_y(" << _x<<")" << std::endl
                << "determinent is negative, " << det << std::endl;
      throw std::runtime_error("");
    }
    det = std::sqrt(det);
    types::t_real u = 1.0 / ( 2.0 * a );  
    types::t_real r0 =  (-b + det ) * u;
    types::t_real r1 =  (-b - det ) * u;
    if ( std::abs(r0 - 1.0 ) < types::tolerance ) r0 = 1.0;
    if ( std::abs(r0 + 1.0 ) < types::tolerance ) r0 = -1.0;
    if ( std::abs(r1 - 1.0 ) < types::tolerance ) r1 = 1.0;
    if ( std::abs(r1 + 1.0 ) < types::tolerance ) r1 = -1.0;
    if (     ( r0 < -1.0 or r0 > 1.0 )
         and ( r1 < -1.0 or r1 > 1.0 ) )
    {
      std::cerr << a + b +c << " " << a - b + c << std::endl;
      std::cerr << "Error when using Concentration::get_y(" << _x<< ")" << std::endl;
      std::cerr << " r0= " << r0  << " and r1= " << r1 << std::endl;
      throw std::runtime_error("");
    }
    if ( r0 < -1.0 or r0 > 1.0 ) 
      return r1; 
    return r0;
  }
//! \deprecated normalizes the concentration of  a structure.
//! Replaced by SingleSite::Concentration and affiliated
types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_real _target = -2.0);


/** @} */
#endif 

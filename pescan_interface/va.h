//
//  Version: $Id$
//
#ifndef _PESCAN_VA_H_
#define _PESCAN_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "interface.h"
#include <vff/functional.h>

#include <opt/types.h>
#include <opt/function_base.h>

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace Pescan
{
  class VirtualAtom : public functional::Base<>
  {
     protected:
       //! Type from which the VA functional is derived
       typedef Functional t_Base;
       //! Type of the pescan interface class
       typedef Functional t_Pescan;
       //! Type of the Valence Force Field Functional class
       typedef Functional t_Vff;

     protected:
       Ising_CE::Structure &structure;
       t_Vff vff;
       t_Pescan pescan;
       Bands result;
       types::t_real deriv_amplitude;


     public:
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_Base(), structure(_str),
                     vff( structure ), pescan(),
                     result(), deriv_amplitude(0.01) {}
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : t_Base( _c ), structure( _c.structure ),
                     vff(_c.structure), pescan( _c.pescan ),
                     result( _c.result), deriv_amplitude( _c.deriv_amplitude ) {}
        

       // Simple constainer behaviors required by Minimizer::VA and
       // Minimizer::Beratan

       //! Returns the size of VirtualAtom::va_vars.
       types::t_unsigned size() const { return va_vars.size(); }
       //! Returns an iterator to the first \e virtual variable (atomic occupation).
       t_Container::iterator begin() { return va_vars.begin(); }
       //! \brief Returns an iterator to one past the last \e virtual variable
       //!        (atomic occupation).
       t_Container::iterator end() { return va_vars.end(); }
       //! \brief Returns a constant iterator to the first \e virtual variable
       //!        (atomic occupation).
       t_Container::const_iterator begin() const { return va_vars.begin(); }
       //! \brief Returns a constant iterator to one past the last \e virtual
       //!        variable (atomic occupation).
       t_Container::const_iterator end() const { return va_vars.end(); }

       // Now truly "functional" stuff.
       
       //! Initializes the variables with respect to Functional::structure.
       bool init();
       //! \brief Evaluated the strain after copying the occupations from
       //!        VirtualAtom::va_vars.
       t_Type evaluate();
       //! Returns the \e virtual gradient in direction \a _pos
       t_Type evaluate_one_gradient( types::t_unsigned _pos );
       //! Computes the \e virtual gradients and returns the energy
       t_Type evaluate_with_gradient( t_Type* _grad );
       //! Computes the \e virtual gradients
       void evaluate_gradient( t_Type* _grad );

       //! Loads pescan and vff minimizers from XML
       bool Load( const TiXmlElement &_node );
         {  return pescan.Load( _node )  and vff.Load( _node ); }

     protected:
       //! Transfers occupations from VirtualAtom::va_vars to Functional::structure.
       void unpack_variables();

       t_Type position_grad( types :: t_int _pos );
       t_Type potential_grad( types :: t_int _pos );

  }


} // namespace Pescan

#endif


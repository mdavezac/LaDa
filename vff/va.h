//
//  Version: $Id$
//
#ifndef _VFF_VA_H_
#define _VFF_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/va_function.h>
#include <lamarck/structure.h>

#ifdef _DOFORTRAN
#include <opt/opt_frprmn.h>
#else
#include <opt/gsl_minimizers.h>
#endif

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace Vff
{

  //! \brief Implements a Virtual Atom functional around Vff::Functional derived object.
  //! \details In other words, this functional is capable of returning the
  //!          gradient with respect to a change in the atomic occupation
  //!          within the structure. It should interact quite well with
  //!          Minimizer::VA and Minimizer::Beratan, Note that the vff
  //!          functional members are hidden.
  //! \note This class is mean to be used with a vff derived object.
  //!       At this two types have been tried:
  //!         - with Vff::Functional as the base class, 
  //!         - with Vff::Layered as the base class.
  //!         .
  //!       You can declare other virtual atom functionals at your own risk.
  template <class T_VFFBASE>
  class VABase : protected T_VFFBASE, public function::VirtualAtom
  {
     //! Type from which the vff functional is derived
     typedef T_VFFBASE t_VffBase;
     //! Type from which the VA functional is derived
     typedef function::VirtualAtom t_VABase;

    protected:
     //! The type of the atom container
     typedef Ising_CE :: Structure :: t_Atoms   t_Atoms;
     //! The type of the atom
     typedef t_Atoms :: value_type t_Atom;
     //! Type of the container holding the atomic centers
     typedef typename t_VffBase :: t_Centers t_Centers;  
     //! Type of the atomic centers
     typedef typename t_Centers :: value_type t_Center;  
     //! Type of the container holding the atomic functionals
     typedef typename t_VffBase :: t_AtomicFunctionals t_AtomicFunctionals;  
     //! Type of the atomic functionals
     typedef typename t_AtomicFunctionals :: value_type t_AtomicFunctional;  


     // Those public members of t_VffBase which are expected in this class
    protected:
     using t_VffBase::energy;

     // Those protected members of t_VffBase which are expected in this class
    protected:
     using t_VffBase::centers;
     using t_VffBase::functionals;


     public:
       //! see functional::Base::t_Type
       typedef t_VABase :: t_Type t_Type;
       //! see functional::Base::t_Container
       typedef t_VABase :: t_Container  t_Container;
       #ifdef _DOFORTRAN
         //! Type of the minimizer for minimizing strain
         typedef Minimizer::Frpr<t_VffBase> t_Minimizer;
       #else
         //! Type of the minimizer for minimizing strain
         typedef Minimizer::GnuSL<t_VffBase> t_Minimizer;
       #endif

     protected:
       //! The minimizer with which vff is minimized
       t_Minimizer minimizer;

     public:
#ifdef _DOFORTRAN
       //! Constructor and Initializer
       VABase   ( Ising_CE::Structure &_str )
              : t_VffBase( _str ), t_VABase( _str ),
                minimizer( *this, vff_for_frprmn )
         { vff_for_frprmn( (double*) this, NULL ); }
#else
       //! Constructor and Initializer
       VABase   ( Ising_CE::Structure &_str )
              : t_VffBase( _str ), t_VABase( _str ), minimizer( *this ) {}
#endif
       //! Copy Constructor
       VABase   ( const VABase &_c )
              : t_VffBase( _c ), t_VABase( _c ), minimizer( *this ) {}
       //! Destructor
       ~VABase() {}
        
       //! Loads the vff's and the minimizer's parameters from XML
       bool Load( const TiXmlElement &_node );

       // Now truly "functional" stuff.
       
       //! \brief Evaluated the strain after copying the occupations from
       //!        VirtualAtom::va_vars.
       t_Type evaluate();
       //! Returns the \e virtual gradient in direction \a _pos
       t_Type evaluate_one_gradient( types::t_unsigned _pos );
       //! Computes the \e virtual gradients and returns the energy
       t_Type evaluate_with_gradient( t_Type* _grad );
       //! Computes the \e virtual gradients
       void evaluate_gradient( t_Type* _grad );
       //! Forwards Vff::Functional::print_escan_input()
       void print_escan_input( const std::string &_f = "atom.config") const
         { t_VffBase::print_escan_input( _f ); }
       //! \brief Initializes the va variables, and optionnally the centers
       //! \details This routine compounds the function::Base::init()
       //!          capabilities expected by Minimizer objects, while adding
       //!          the capacity for recomputing the first-neighbor tree of
       //!          vff.
       bool init( bool _redocenters = false )
        { return  t_VABase::init() and ( _redocenters ? t_VffBase::construct_centers(): true ); }

     protected:

  };

} // namespace vff 

#include "va.impl.h"

#endif // _VFF_FUNCTIONAL_H_

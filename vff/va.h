//
//  Version: $Id$
//
#ifndef _LADA_VFF_VA_H_
#define _LADA_VFF_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/path.hpp>
#include <boost/mpl/vector.hpp>

#include <opt/va_function.h>
#include <crystal/structure.h>
#include <mpi/mpi_object.h>

#include "functional.h"

#include <minimizer/function_wrapper.h>
#include <minimizer/any.h>
#include <minimizer/frprmn.h>
#include <minimizer/gsl_mins.h>
#include <minimizer/decoupled.h>


namespace LaDa
{
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
       //! Type of the path.
       typedef boost::filesystem::path t_Path;
      public:
       //! Type from which the vff functional is derived
       typedef T_VFFBASE t_VffBase;
       //! Type from which the VA functional is derived
       typedef function::VirtualAtom t_VABase;

      protected:
       //! The type of the atom container
       typedef Crystal :: Structure :: t_Atoms   t_Atoms;
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
         //! Type of the minimizer for minimizing strain
         typedef Minimizer::Any
                 < 
                   boost::mpl::vector
                   <
                     Minimizer::Frpr, 
                     Minimizer::Gsl, 
                     Minimizer::Decoupled
                   > 
                 > t_Minimizer;

       protected:
         //! The minimizer with which vff is minimized
         t_Minimizer minimizer;

       public:
         //! Constructor and Initializer
         VABase   ( Crystal::Structure &_str )
                : t_VffBase( _str ), t_VABase( _str ), minimizer() {}
         //! Copy Constructor
         VABase   ( const VABase &_c )
                : t_VffBase( _c ), t_VABase( _c ), minimizer( _c.minimizer ) {}
         //! Destructor
         ~VABase() {}
          
         //! Loads the vff's and the minimizer's parameters from XML
         bool Load( const TiXmlElement &_node );

         // Now truly "functional" stuff.
         
         //! \brief Evaluated the strain after copying the occupations from
         //!        VirtualAtom::va_vars.
         t_Type evaluate();
      //  //! Returns the \e virtual gradient in direction \a _pos
      //  t_Type evaluate_one_gradient( types::t_unsigned _pos );
      //  //! Computes the \e virtual gradients and returns the energy
      //  t_Type evaluate_with_gradient( t_Type* _grad );
      //  //! Computes the \e virtual gradients
      //  void evaluate_gradient( t_Type* _grad );
         //! Forwards Vff::Functional::print_escan_input()
         void print_escan_input( const t_Path& _f = "atom.config") const
           { t_VffBase::print_escan_input( _f ); }
         //! \brief Initializes the va variables, and optionnally the centers
         //! \details This routine compounds the function::Base::init()
         //!          capabilities expected by Minimizer objects, while adding
         //!          the capacity for recomputing the first-neighbor tree of
         //!          vff.
         bool init( bool _redocenters = false );
         //! Returns reference to the Vff base class
         t_VffBase& Vff() { return *static_cast<t_VffBase*>(this); }
         //! Returns constant reference to the Vff base class
         const t_VffBase& Vff() const
          { return *static_cast<const t_VffBase*>(this); }
         //! gets already computed stress from vff. 
         void get_stress( atat::rMatrix3d &_s ) const { _s = Vff().stress; }

       protected:

    };

  } // namespace vff 
} // namespace LaDa
#include "va.impl.h"

#endif // _VFF_FUNCTIONAL_H_

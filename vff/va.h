#ifndef _LADA_VFF_VA_H_
#define _LADA_VFF_VA_H_

#include "LaDaConfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/mpl/vector.hpp>

#include <opt/va_function.h>
#include <crystal/structure.h>
#include <mpi/mpi_object.h>

#include "functional.h"

#include <minimizer/function_wrapper.h>
#include <minimizer/frprmn.h>
#ifdef LADA_WITH_GSL
# include <minimizer/gsl_mins.h>
#endif
#ifdef LADA_WITH_MINUIT2
# include <minimizer/minuit2.h>
#endif
#if defined(LADA_WITH_GSL) or defined(LADA_WITH_MINUIT2)
# include <minimizer/variant.h>
#endif
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
        
#       if defined(LADA_WITH_GSL) or defined(LADA_WITH_MINUIT2)
        //! Type of the minimizer for minimizing strain
        typedef Minimizer::Variant
                < 
                  boost::mpl::vector
                  <
                    Minimizer::Frpr
#                   ifdef LADA_WITH_GSL
                      , Minimizer::Gsl 
#                   endif
#                   ifdef LADA_WITH_MINUIT2
                      , Minimizer::Minuit2
#                   endif
                  > 
                > t_Minimizer;
#       else
          //! Type of the minimizer for minimizing strain
          typedef Minimizer::Frpr t_Minimizer;
#       endif

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
        t_Type evaluate(bool relax=true);
        //! Returns the \e virtual gradient in direction \a _pos
        t_Type evaluate_one_gradient( types::t_unsigned _pos );
        //! Computes the \e virtual gradients and returns the energy
        t_Type evaluate_with_gradient( t_Type* _grad );
        //! Computes the \e virtual gradients
        void evaluate_gradient( t_Type* _grad );
        //! Forwards Vff::Functional::print_escan_input()
        void print_escan_input( const t_Path& _f = "atom.config") const
          { t_VffBase::print_escan_input( _f ); }
        //! \brief Initializes the va variables, and optionnally the centers
        //! \details This routine compounds the function::Base::init()
        //!          capabilities expected by Minimizer objects, while adding
        //!          the capacity for recomputing the first-neighbor tree of
        //!          vff.
        bool init( bool _redocenters = false, bool _verbose = false );
        //! Returns reference to the Vff base class
        t_VffBase& Vff() { return *static_cast<t_VffBase*>(this); }
        //! Returns constant reference to the Vff base class
        const t_VffBase& Vff() const
         { return *static_cast<const t_VffBase*>(this); }
        //! Returns constant structure.
        const Crystal::Structure& get_structure() const { return structure; }
        

         //! \brief Returns a reference to the computed stress
         //! \sa Functional::stress
         const math::rMatrix3d& get_stress() const { return t_VffBase::stress; }

#        ifdef LADA_MPI
           //! Sets mpi pointer.
           void set_mpi( boost::mpi::communicator* _c )
             { t_VffBase ::  set_mpi( _c ); }
           //! Returns reference to communicator.
           boost::mpi::communicator &comm() { return t_VffBase :: comm(); }
           //! Returns a constant reference to communicator.
           const boost::mpi::communicator &comm() const { return t_VffBase :: comm(); } 
#        endif

         //! Sets the bond parameters.
         template< class T_TUPLE >
           void set_bond( const std::string &_type, const T_TUPLE& _tuple )
             { t_VffBase::set_bond( _type, _tuple ); }
         //! Returns bond parameters, first the length, then the alphas.
         boost::tuples::tuple< const types::t_real&, const types::t_real&,
                               const types::t_real&, const types::t_real&, 
                               const types::t_real&, const types::t_real& >
           get_bond( const std::string &_type ) const
             { return t_VffBase::get_bond( _type ); }
         //! Sets the angle parameters.
         template< class T_TUPLE >
           void set_angle( const std::string &_type, const T_TUPLE& _tuple )
             { t_VffBase::set_angle( _type, _tuple); }
         //! Returns angle parameters, first the length, then sigma, then the betas.
         boost::tuples::tuple< const types::t_real&, const types::t_real&, 
                               const types::t_real&, const types::t_real&,
                               const types::t_real&, const types::t_real&,
                               const types::t_real& >
           get_angle( const std::string &_type ) const
             { return t_VffBase::get_angle( _type ); }

         //! Copies parameters from argument.
         template<class T> 
           void copy_parameters(VABase<T> const &_f) { t_VffBase::copy_parameters(_f.Vff()); }

         //! Sets minimizer. For python mostly.
         void set_minimizer(t_Minimizer const &_minimizer) { minimizer = _minimizer; }

      protected:
        //! The minimizer with which vff is minimized
        t_Minimizer minimizer;
    };

  } // namespace vff 
} // namespace LaDa
#include "va.impl.h"

#endif // _VFF_FUNCTIONAL_H_

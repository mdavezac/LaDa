#ifndef _LADA_VFF_VA_H_
#define _LADA_VFF_VA_H_

#include "LaDaConfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/mpl/vector.hpp>

#include <opt/va_function.h>
#include <crystal/structure.h>
#include <opt/mpi.h>

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
  namespace vff
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
    template <class T_VFFBASE>
    class WithMinimizer : protected T_VFFBASE
    {
      //! Type of the path.
      typedef boost::filesystem::path t_Path;
      public:
        //! Type from which the vff functional is derived
        typedef T_VFFBASE t_VffBase;
        
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
        //! The minimizer with which vff is minimized
        t_Minimizer minimizer;

        //! Type of the return.
        typedef typename t_VffBase::t_Return t_Return;
      protected:
        //! The type of the atom container
        typedef typename Vff::t_Arg::t_Atoms   t_Atoms;
        //! The type of the atom
        typedef t_Atoms :: value_type t_Atom;
        //! Type of the container holding the atomic centers
        typedef typename t_VffBase :: t_Centers t_Centers;  
        //! Type of the atomic centers
        typedef typename t_Centers :: value_type t_Center;  


        // Those public members of t_VffBase which are expected in this class
      public:
        using t_VffBase::check_input;
      protected:
        using t_VffBase::energy;

        // Those protected members of t_VffBase which are expected in this class
      protected:
        using t_VffBase::centers_;
        using t_VffBase::structure;

      public:
        //! Constructor and Initializer
        WithMinimizer() {}
        //! Constructor and Initializer
        WithMinimizer   ( typename T_VFFBASE::t_Arg &_str )
               : t_VffBase( _str ), minimizer() {}
        //! Copy Constructor
        WithMinimizer   ( const WithMinimizer &_c )
               : t_VffBase( _c ), minimizer( _c.minimizer ) {}
        //! Destructor
        ~WithMinimizer() {}
         
        // Now truly "functional" stuff.
        
#       ifdef LADA_MPI
          //! \brief Evaluated the strain after copying the occupations from
          //!        VirtualAtom::va_vars.
          t_Return evaluate(boost::mpi::communicator const &_comm, bool relax=true)
            { T_VFFBASE::comm = _comm; return evaluate(relax); }
      protected:
#       endif
        //! \brief Evaluated the strain after copying the occupations from
        //!        VirtualAtom::va_vars.
        t_Return evaluate(bool relax=true);
#       ifdef LADA_MPI
      public:
#       endif
//       //! Returns the \e virtual gradient in direction \a _pos
//       t_Type evaluate_one_gradient( types::t_unsigned _pos );
//       //! Computes the \e virtual gradients and returns the energy
//       t_Type evaluate_with_gradient( t_Type* _grad );
//       //! Computes the \e virtual gradients
//       void evaluate_gradient( t_Type* _grad );
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
        t_VffBase& VffBase() { return *static_cast<t_VffBase*>(this); }
        //! Returns constant reference to the Vff base class
        const t_VffBase& VffBase() const
         { return *static_cast<const t_VffBase*>(this); }
        //! Returns constant structure.
        Crystal::TStructure<std::string> const & get_structure() const { return structure; }
        //! Sets internal structure.
        void set_structure(Crystal::TStructure<std::string> const & _str) { structure = _str; }
        

         //! \brief Returns a reference to the computed stress
         //! \sa Functional::stress
         const math::rMatrix3d& get_stress() const { return t_VffBase::stress; }

         //! Sets the bond parameters.
         template< class T_TUPLE >
           void set_bond( const std::string &_type, const T_TUPLE& _tuple )
             { t_VffBase::set_bond( _type, _tuple ); }
         //! Copies parameters from argument.
         template<class T> 
           void copy_parameters(WithMinimizer<T> const &_f) { t_VffBase::copy_parameters(_f.VffBase()); }
    };

    template< class T_VFFBASE > typename WithMinimizer<T_VFFBASE> :: t_Return 
      WithMinimizer<T_VFFBASE> :: evaluate(bool relax)
      {
        // no minimization required if variables is empty.
        typename t_VffBase :: t_Arg arg;
        t_VffBase :: init( arg );
          
        if(arg.size() and relax) minimizer( *( (t_VffBase*) this), arg );
     
        t_VffBase :: structure.energy = t_VffBase::energy();

        return t_VffBase::structure.energy;
      }

    template< class T_VFFBASE > 
      inline bool WithMinimizer<T_VFFBASE> :: init(bool _redocenters, bool _verbose)
      {
        if( not _redocenters) return true;
        return t_VffBase::initialize_centers(_verbose);
      }

  } // namespace vff 
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_

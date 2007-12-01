//
//  Version: $Id$
//
#ifndef _VFF_VA_H_
#define _VFF_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "functional.h"

#include <opt/va_function.h>

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

  //! \brief Implements a Virtual Atom functional around Vff::Functional.
  //! \details In other words, this functional is capable of returning the
  //!          gradient with respect to a change in the atomic occupation
  //!          within the structure. It should interact quite well with
  //!          Minimizer::VA and Minimizer::Beratan, Note that the vff
  //!          functional members are hidden.
  class VirtualAtom : protected Functional, public function::VirtualAtom
  {
     //! Type from which the vff functional is derived
     typedef Functional t_VffBase;
     //! Type from which the VA functional is derived
     typedef function::VirtualAtom t_VABase;

    protected:
     //! The type of the atom container
     typedef Ising_CE :: Structure :: t_Atoms   t_Atoms;
     //! The type of the atom
     typedef t_Atoms :: value_type t_Atom;
     //! Type of the container holding the atomic centers
     typedef t_VffBase::t_Centers t_Centers;  
     //! Type of the atomic centers
     typedef t_Centers::value_type t_Center;  
     //! Type of the container holding the atomic functionals
     typedef t_VffBase::t_AtomicFunctionals t_AtomicFunctionals;  
     //! Type of the atomic functionals
     typedef t_AtomicFunctionals::value_type t_AtomicFunctional;  



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
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_VffBase( _str ), t_VABase( _str ), minimizer( *this, vff_for_frprmn )d
         { vff_for_frprmn( (double*) this, NULL ); }
#else
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_VffBase( _str ), t_VABase( _str ), minimizer( *this ) {}
#endif
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : t_VffBase( _c ), t_VABase( _c ), minimizer( *this ) {}
       //! Destructor
       ~VirtualAtom() {}
        
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

     protected:

  };

  /** \brief Adds a positional derivate for Virtual Atom with pescan. 
   *  \details More specifically, the escan Hamiltonian can be derived as (see
   *           Pesca::VirtualAtom),
   *           \f[
   *             \frac{\mathrm{d}\mathcal{H}}{\mathrm{d}\mathcal{S}} = 
   *               \underbrace{
   *                 \sum_{\mathbf{R}} \frac{\partial\mathcal{H}}
   *                                        {\partial\alpha_{\mathbf{R}}}
                                       \frac{\partial\alpha_{\mathbf{R}}}
   *                                        {\partial S}
   *               }_{\Delta \mathcal{H}_{\mathrm{chem}}}
   *               +
   *               \underbrace{
   *                 \sum_{\mathbf{R}}\nabla_{\mathbf{R}}\mathcal{H}\cdot
   *                   \frac{\partial\mathcal{R}}{\partial S}
   *               +
   *                 \sum_{\mathbf{R}}\frac{\partial \mathcal{H} }
   *                                       {\partial\epsilon_\mathcal{R}}
   *                                  \frac{\partial\epsilon_\mathcal{R}}
   *                                       {\partial S}
   *               }_{\Delta \mathcal{H}_{\mathrm{stress}}}
   *               +
   *               \underbrace{
   *                 \frac{\partial \mathcal{H}}{\partial\Omega}
   *                 \frac{\partial \Omega}{\partial S}
   *               }_{\Delta \mathcal{H}_{\mathrm{cell}}}
   *           \f]
   *           This class is capable of writing the atomic input of both
   *           perturbed potentials \f$\Delta \mathcal{H}_{\mathrm{chem}}\f$
   *           and \f$\Delta \mathcal{H}_{\mathrm{stress}}\f$. It is called
   *           upon by Pescan::VirtualAtom.
   */
  class PescanPerturbations : public VirtualAtom
  {
     //! Helps in determining escan pseudo
     typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
     //! Helps in determining escan pseudo
     typedef std::vector< t_pseudo > t_pseudos;
     //! The type of the atom container
     typedef Ising_CE::Structure::t_Atoms t_Atoms;
     //! The type of the atom
     typedef Ising_CE::Structure::t_Atom  t_Atom;

    public:
     //! The amplitude of the numerical derivative
     const static types::t_real deriv_amplitude = 0.01;

    protected:
      //! Type of the container holding the atomic centers
      typedef VirtualAtom::t_Centers t_Centers;  
      //! Type of the atomic centers
      typedef t_Centers::value_type t_Center;  
      //! Type of the container holding the atomic functionals
      typedef VirtualAtom::t_AtomicFunctionals t_AtomicFunctionals;  
      //! Type of the atomic functionals
      typedef t_AtomicFunctionals::value_type t_AtomicFunctional;  
      //! Type of the VA functional
      typedef function :: VirtualAtom t_VA;

    public:
      //! Base class for this class. 
      typedef VirtualAtom t_Base;

    public:
      //! File to which to print atomic configurations
      std::string filename;

    public:
      //! Constructor and Initializer
      PescanPerturbations   ( Ising_CE::Structure &_str )
                         : t_Base( _str ), filename("atomic_config") {}
      //! Copy Constructor
      PescanPerturbations   ( const PescanPerturbations &_c )
                         : t_Base( _c ), filename(_c.filename) {}
      //! Destructor
      ~PescanPerturbations() {}
    
      //! Printes \f$\Delta \mathcal{H}_{\mathrm{chem}}\f$ input for Pescan::VirtualAtom
      void chemical( types::t_unsigned _pos );
      //! Printes \f$\Delta \mathcal{H}_{\mathrm{stress}}\f$ input for Pescan::VirtualAtom
      types::t_real stress( types::t_unsigned _pos );
      //! Printes \f$\mathcal{H}\f$ input to file for PescanVirtualAtom
      void zero_order() { t_Base :: print_escan_input( filename ); }

    protected:
      //! \brief Finds escan potentials
      //! \details Escan potentials correspond to a mixing of the potential of
      //!          A interacting with each of its first neighbors. This
      //!          function finds how many types there are of first neighbors.
      void find_escan_pseudos( const t_Center &_center,
                               t_pseudos _pseudos );
  };

} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_

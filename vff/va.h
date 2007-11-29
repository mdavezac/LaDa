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
#include <opt/gsl_minimizers.h>

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace Vff
{

  //! \brief Implements a Virtual Atom functional around Vff::Functional.
  //! \details In other words, this functional is capable of returning the
  //!          gradient with respect to a change in the atomic occupation
  //!          within the structure. It should interact quite well with
  //!          Minimizer::VA and Minimizer::Beratan,
  class VirtualAtom : protected Functional, public function::VirtualAtom
  {
     //! Type from which the VA functional is derived
     typedef Functional t_VffBase;
     //! Type from which the VA functional is derived
     typedef Functional t_VABase;

    protected:
     //! The type of the atom container
     typedef t_VABase :: t_Atoms   t_Atoms;
     //! The type of the atom
     typedef t_Atoms :: value_type t_Atom;


     public:
       //! see functional::Base::t_Type
       typedef t_VABase :: t_Type t_Type;
       //! see functional::Base::t_Container
       typedef t_VABase :: t_Container  t_Container;
       //! Type of the minimizer for minimizing strain
       typedef Minimizer::GnuSL<t_VffBase> t_Minimizer;

     protected:
       //! The minimizer with which vff is minimized
       t_Minimizer minimizer;

     public:
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_VffBase( _str ), t_VABase( _str ), minimizer( *this ) {}
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : t_VffBase( _c ), t_VABase( _c ), minimizer( *this ) {}
        
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

     protected:
  };

  /** \brief Adds a positional derivate for Virtual Atom with pescan. 
   *  \details More specifically, any Virtual Atom gradient of a physical
   *           property will include a partial derivative with respect to the
   *           configuration and with respect to the position. More
   *           specifically, if \f$\mathcal{O}(\sigma)\f$ is the property,
   *           \f$\sigma=\{S_0,\dots,S_N\}\f$ the configuration, and
   *           \f$\{\hat{R}\}_{\mathrm{min}}\f$ the \e relaxed postions of the
   *           atoms, then:
   *           \f[
   *               \frac{\mathrm{d}\mathcal{O}(\sigma)}{\mathrm{d}S_i} = 
   *                 \left.\frac{\partial \mathcal{O}(\sigma)}{\partial
   *                 S_i}\right|_{\{\hat{R}_{\mathrm{min}}\}}
   *                 + \frac{\partial \mathcal{O}(\sigma)}{\partial \{\hat{R}_{min}\}}
   *                 \cdot \frac{\partial \{\hat{R}_{min}\}}{\partial S_i}
   *           \f]
   *           This partical class takes care of computing the variation of the
   *           positions for a certain change in the occupation (see
   *           PescanPosGrad::deriv_amplitude). It also writes an file for
   *           inputing into escan the resulting \e perturbation. 
   *
   *           The variation of the positions is computed by minimizing the
   *           structure with a flipped atom.
   */
  class PescanPosGrad : public VirtualAtom
  {
     //! Helps in determining escan pseudo
     typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
     //! Helps in determining escan pseudo
     typedef std::vector< t_pseudo > t_pseudos;
     //! The type of the atom container
     typedef Ising_CE::Structure::t_Atoms t_Atoms;
     //! The type of the atom
     typedef Ising_CE::Structure::t_Atom  t_Atom;

     //! The amplitude of the numerical derivative
     const static types::t_real deriv_amplitude = 0.01;

    protected:
      //! Type of the container holding the atomic centers
      typedef typename t_Base::t_Centers t_Centers;  
      //! Type of the atomic centers
      typedef typename t_Centers::value_type t_Center;  
      //! Type of the container holding the atomic functionals
      typedef typename t_Base::t_AtomicFunctionals t_AtomicFunctionals;  
      //! Type of the atomic functionals
      typedef typename t_AtomicFunctionals::value_type t_AtomicFunctional;  


    public:
      //! Base class for this class. 
      typedef VirtualAtom t_Base;

    public:
      //! Constructor and Initializer
      PescanPosGrad( Ising_CE::Structure &_str ) : t_Base( &_str ) {}
      //! Constructor and Initializer
      PescanPosGrad( const t_Base &_c ) : t_Base( &_c ) {}
      //! Copy Constructor
      PescanPosGrad( const PescanPosGrad &_c ) : t_Base( &_c ) {}
    
      //! \brief returns the spin flip energy for spin \a _pos.
      //! \details Also writes to file a pertubed pescan potential
      //!          corresponding to PescanPosGrad::deriv_amplitude\% variation
      //!          of the positions from the starting position to the structure
      //!          with a flipped atom.
      types::t_real position_gradient( types::t_unsigned _pos );
      types::t_real potential_gradient( types::t_unsigned _pos );

    protected:
      //! \brief Finds escan potentials
      //! \details Escan potentials correspond to a mixing of the potential of
      //!          A interacting with each of its first neighbors. This
      //!          function finds how many types there are of first neighbors.
      void find_escan_pseudos( typename t_Centers::const_iterator &_i_center,
                               t_pseudos _pseudos );
  };

} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_

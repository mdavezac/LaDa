#ifndef _VFF_PESCAN_PERTURBATION_H_
#define _VFF_PESCAN_PERTURBATION_H_

#include "LaDaConfig.h"


#include <opt/va_function.h>
#include <boost/filesystem/path.hpp>

#include "functional.h"
#include "va.h"

#include <mpi/mpi_object.h>

namespace LaDa
{
  namespace Vff
  {
    /** \brief Adds a positional derivate for Virtual Atom with pescan. 
     *  \details More specifically, the escan Hamiltonian can be derived as (see
     *           Pescan::VirtualAtom),
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
    template< class T_VFFBASE >
    class PescanPerturbations : public VABase<T_VFFBASE>
    {
        //! Type of the path.
        typedef boost::filesystem::path t_Path;
      public:
        //! Vff base class.
        typedef T_VFFBASE t_Functional;
        //! Base class for this class. 
        typedef VABase<t_Functional> t_Base;
        //! The amplitude of the numerical derivative
        const static types::t_real deriv_amplitude = 0.01;

      protected:
       //! Helps in determining escan pseudo
       typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
       //! Helps in determining escan pseudo
       typedef std::vector< t_pseudo > t_pseudos;
       //! The type of the atom container
       typedef Crystal::Structure::t_Atoms t_Atoms;
       //! The type of the atom
       typedef Crystal::Structure::t_Atom  t_Atom;

       //! Type of the container holding the atomic centers
       typedef typename t_Base :: t_Centers t_Centers;  
       //! Type of the atomic centers
       typedef typename t_Base :: t_Center t_Center;  
       //! Type of the container holding the atomic functionals
       typedef typename t_Base :: t_AtomicFunctionals t_AtomicFunctionals;  
       //! Type of the atomic functionals
       typedef typename t_Base :: t_AtomicFunctional t_AtomicFunctional;  
       //! Type of the VA functional
       typedef typename t_Base :: t_VABase t_VA;


      public:
        //! File to which to print atomic configurations
        t_Path filename;

      public:
        //! Constructor and Initializer
        PescanPerturbations   ( Crystal::Structure &_str )
                           : t_Base( _str ), filename("atomic_config") {}
        //! Copy Constructor
        PescanPerturbations   ( const PescanPerturbations &_c )
                           : t_Base( _c ), filename(_c.filename) {}
        //! Destructor
        ~PescanPerturbations() {}
      
        //! Prints \f$\Delta \mathcal{H}_{\mathrm{chem}}\f$ input for Pescan::VirtualAtom
        void chemical( types::t_unsigned _pos );
        //! Prints \f$\Delta \mathcal{H}_{\mathrm{stress}}\f$ input for Pescan::VirtualAtom
        types::t_real stress( types::t_unsigned _pos );
        //! Prints \f$\mathcal{H}\f$ input to file for PescanVirtualAtom
        void zero_order( const t_Path &_f )
          { t_Base :: print_escan_input( _f ); }
       
        //! gets already computed stress from vff. 
        void get_stress( math::rMatrix3d &_s ) const { _s = t_Base::Vff().get_stress(); }
     
        using t_Base :: evaluate;

      protected:
        //! \brief Finds escan potentials
        //! \details Escan potentials correspond to a mixing of the potential of
        //!          A interacting with each of its first neighbors. This
        //!          function finds how many types there are of first neighbors.
        void find_escan_pseudos( const t_Center &_center,
                                 t_pseudos _pseudos );

        using t_Base :: centers;
        using t_Base :: functionals;
        using t_Base :: structure0;
    };


  } // namespace vff 

} // namespace LaDa
#include "pescan_perturbation.impl.h"

#endif // _VFF_PESCAN_PERTURBATION_H_

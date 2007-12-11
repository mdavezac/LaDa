//
//  Version: $Id$
//
#ifndef _PESCAN_VA_H_
#define _PESCAN_VA_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bandgap.h"

#include <vff/pescan_perturbation.h>
#include <opt/va_function.h>

#include <opt/types.h>

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace Pescan
{
  /** \brief Implements Virtual Atom for escan band-gaps.
   *  \details The Hamiltonian of escan in reciprocal space, with \b G and \b G' two reciprocal vectors is
   *           \f[
   *              \mathcal{H}_{\mathbf{G}, \mathbf{G}'} =
   *              \frac{1}{2}G^2\delta_{\mathbf{G}, \mathbf{G}'} +
   *              \mathcal{V}(\mathbf{G}-\mathbf{G}')
   *           \f] 
   *           where the first term is the kinetic energy, and the second the
   *           Fourier transform of hte screened pseudo-potentials.
   *           The screened-potential is a sum over the
   *           (real-space) lattice sites \b R of spherical potentials and of
   *           strain-dependent terms,
   *           \f[
   *               \mathcal{V}(\mathbf{r}) = \sum_{\mathbf{R}, \alpha_{\mathbf{R}} } 
   *                 v_{\alpha_{\mathbf{R}}}( |\mathbf{r}-\mathbf{R}| ) [ 1 +
   *                 \gamma_{\alpha_{\mathbf{R}}} \mathrm{Tr}(\epsilon) ]
   *           \f]
   *           with \f$\alpha_{\mathbf{R}}\f$ denoting the dependence with respect to the
   *           chemical specie at site \b R and its immediate surrounding,
   *           \f$v_{\alpha_{\mathbf{R}}}\f$ the spherical atomic potential,
   *           \f$\gamma_{\alpha_{\mathbf{R}}}\f$ the dependence with respect to
   *           the microscopic strain, and \f$\mathrm{Tr}(\epsilon)\f$ the
   *           trace of the strain tensor \f$\epsilon\f$ (the dependence with
   *           respect to the microscopic strain is expressed for zing-blende
   *           structures here).  Finally, the dependence with respect to the
   *           local environnment is as follows, for \f$A_0\f$ atom surrounded
   *           by a tetrahedron of \f$p B_0\f$ atoms and \f$(4-p) B_1\f$ atoms,
   *           \f[
   *             v_{A_0}[ p B_0,(4-p) B_1] = \frac{4-p}{4} v_{A_0}[A_0B_1]
   *             \frac{p}{4} v_{A_0}[A_0B_0] 
   *           \f]
   *           In all the above, the \e dependence of the \b position \b
   *           vectors \b R  and the \e microscopic \e strains \f$\epsilon\f$
   *           \e with \e respect \e to \e the chemical \e species \e is \e
   *           implicit. Indeed, each configuration must be
   *           strain-relaxed with (generally) the Valence Force Field method.
   *
   *           Eventually, the Hamiltonian can be expressed as
   *           \f$\mathcal{H}\left[\sigma,
   *           \{\mathbf{R}\}_{\mathrm{min}},
   *           \epsilon_{\mathbf{R}}, \Omega\right]\f$, with
   *           \f$\sigma\f$ the configuration, 
   *           \f$\{\mathbf{R}\}_{\mathrm{min}}\f$ the
   *           strain-relaxed fractional atomic positions, 
   *           \f$\epsilon_{\mathbf{R}}\f$ the relaxed
   *           microscopic strain at \b R, and \f$\Omega\f$ the unit cell. It
   *           is paramount to note that both positions and microscopic strains
   *           are \e input to pescan (generally from Vff). As such, a complete
   *           derivative with respect to a lattice occupation \f$S\f$ at some
   *           site which would include the dependence with respect to the two
   *           geometric factors takes the following form
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
   *           These four terms can be approximated in the way described below.
   *           Let \f$\sigma\f$ be the original configuration, and
   *           \f$\sigma'\f$ the configuration with the \e flipped
   *           spin/occupation at one lattice site. At the start of the
   *           Virtual-Atom gradient occupation, we have:  
   *             - Unit-cell of \f$\sigma\f$, strain-relaxed
   *             - Fractional coordinates of \f$\sigma\f$, strain-relaxed
   *             - Microscopic-strain of \f$\sigma\f$, strain-relaxed
   *             - Wavefunction of \f$\sigma\f$
   *             - Eigenvalues/Bandgap of \f$\sigma\f$
   *             .
   *           The following quantities can be easily obtained from a
   *           (relatively fast) Vff computation and can be used:
   *             - Unit-cell of \f$\sigma'\f$, strain-relaxed
   *             - Fractional coordinates of \f$\sigma'\f$, strain-relaxed
   *             - Microscopic-strain of \f$\sigma'\f$, strain-relaxed
   *             .
   *           The following quantities should be compuationally expensive and
   *           should be avoided:
   *             - Wavefunction of \f$\sigma'\f$
   *             - Eigenvalues/Bandgap of \f$\sigma'\f$
   *             .
   *            Eventually, we can adopt the following basic strategy:
   *             - Compute the pertubation potential induced in one variable
   *             when flipping one spin and keeping all other variables
   *             constant. 
   *             - Apply the unperturbed wavefunctions to obtain the first order eigenvalues.
   *             .
   *           <DL>
   *             <DT> \f$\Delta \mathcal{H}_{\mathrm{chem}}\f$ </DT>
   *             <DD> 
   *                   - \b Keep the wavefunctions, positions, microscopic
   *                     strain, and unit-cell from \f$\sigma\f$
   *                   - \b Construct the perturbed potential. Only five
   *                   lattice-sites are concerned: the lattice-site where the
   *                   occupation is changing and its surrounding tetrahedron.
   *                   The perturbed potential consists of the potential of \f$\sigma\f$
   *                    minus a small part of the potential of \f$\sigma'\f$
   *                   - \b Apply the wave-functions \e as are onto the
   *                   perturbed potential and obtain new band-gap.  Eg escan
   *                   with 0 iterations.
   *                   - See Vff::PescanPertubation::chemical()
   *                   .
   *             </DD>
   *             <DT> \f$\Delta \mathcal{H}_{\mathrm{stress}}\f$
   *             <DD> 
   *                   - \b Keep the wavefunctions, and unit-cell from \f$sigma\f$
   *                   - \b Compute the strain relaxed Vff positions for
   *                     \f$\sigma'\f$ with constant unit-cell and the new microstrain.
   *                   - \b Construct a perturbed potential consisting of the
   *                   old positions/old microstrain - displaced atoms/new microstrain.
   *                   - \b Apply the wave-functions \e as are onto the
   *                   perturbed potential.  Eg escan with 0 iterations.
   *                   - See Vff::PescanPertubation::stress()
   *                   .
   *             </DD>
   *             <DT> \f$\Delta \mathcal{H}_{\mathrm{cell}}\f$
   *             <DD> 
   *                   Hopefully. this term is small... Indeed, once the
   *                   unit-cell is changed, nothing computed from \f$\sigma\f$ can be
   *                   used since the meshes in real-space and reciprocal-space
   *                   of the Fourier transform used in pescan will not match
   *                   between the old unit-cell and the next. 
   *             </DD>
   *           </DL>
   *
   *
   */
  class VirtualAtom : protected BandGap, public function::VirtualAtom
  {
     friend std::ostream& operator<<( std::ostream& _stream, const VirtualAtom& _va );
     protected:
       //! Type from which the VA functional is derived
       typedef function::VirtualAtom t_VABase;
       //! Type of the pescan interface class
       typedef BandGap t_PescanBase;
       //! \brief Type of the Valence Force Field Functional class
       //! \details It is fitted into a class which adds positional pertubation
       //!          stuff
       typedef Vff::PescanPerturbations t_Vff;

       //! Amplitude of the numerical derivative
       const static types::t_real deriv_amplitude = 0.01;

       //! The type of atom container
       typedef Ising_CE::Structure :: t_Atoms t_Atoms;
       //! The type of atoms 
       typedef t_Atoms :: value_type  t_Atom;
       
     public:
       //! For doing chemical gradients
       const static types::t_unsigned CHEMICAL_GRADIENT          = 1;
       //! For doing stress gradients                     
       const static types::t_unsigned STRESS_GRADIENT            = 2;
       //! For doing both chemical and stress gradients
       const static types::t_unsigned CHEMICAL_STRESS_GRADIENTS  = 3;

     protected:
       //! the vff object to get the variation of the positions
       t_Vff vff;
       //! Tracks which gradients to do
       types::t_unsigned do_gradients;
       


     public:
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_PescanBase(), t_VABase( _str ),
                     vff( structure ), do_gradients(CHEMICAL_STRESS_GRADIENTS) 
       {
#ifdef _MPI
         std::ostringstream sstr;
         sstr << vff.atomic_config << "." << mpi::main.rank();
         vff.atomic_config = sstr.str();
#endif
       }
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : t_PescanBase( _c ), t_VABase( _c ),
                     vff(_c.vff), do_gradients(_c.do_gradients) {}
        

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
       bool Load( const TiXmlElement &_node )
         {  return t_PescanBase::Load( _node )  and vff.Load( _node ); }

       //! Returns a reference to the virtual atom vff minimizer
       t_Vff& Vff() { return vff; }
       //! Returns a constant reference to the virtual atom vff minimizer
       const t_Vff& Vff() const { return vff; }
       //! Returns a reference to the BandGap base
       operator BandGap&() { return *( (BandGap*) this ); }
       //! Returns a constant reference to the BandGap base
       operator const BandGap&() const { return *( (const BandGap*) this ); }

     protected:
         types::t_real apply_wfns();
  };

  inline VirtualAtom::t_Type VirtualAtom::evaluate()
  { 
    vff.evaluate();
    vff.zero_order();
    t_PescanBase::do_input_wavefunctions = false;
    t_PescanBase::escan.wavefunction_out = "zero_order";
    
    return t_PescanBase::operator()( structure ); 
  }

  inline VirtualAtom::t_Type VirtualAtom::evaluate_with_gradient( t_Type* _grad )
  { 
    types :: t_real result = t_PescanBase::operator()( structure ); 
    evaluate_gradient( _grad );
    return result;
  }

  inline void VirtualAtom::evaluate_gradient( t_Type* _grad )
  {
    t_Container :: iterator i_var = t_VABase::va_vars.begin();
    t_Container :: iterator i_var_end = t_VABase::va_vars.end();
    t_Type* i_grad = _grad;
    for(types::t_unsigned n=0; i_var != i_var_end; ++i_var, ++i_grad, ++n )
      *i_grad += evaluate_one_gradient( n );
  }

  inline VirtualAtom::t_Type VirtualAtom::apply_wfns()
  { 
    vff.evaluate();
    t_PescanBase::do_input_wavefunctions = true;
    t_PescanBase::escan.wavefunction_in = "zero_order";
    t_PescanBase::escan.wavefunction_out = "first_order";

    types::t_int niter      = t_PescanBase::escan.niter;
    types::t_int nlines     = t_PescanBase::escan.nlines;
    types::t_real tolerance = t_PescanBase::escan.tolerance;

    t_PescanBase::escan.niter     = 0;
    t_PescanBase::escan.nlines    = 0;        
    t_PescanBase::escan.tolerance = 10000;    
    
    types::t_real result = t_PescanBase::operator()( structure ); 

    t_PescanBase::escan.niter     = niter;
    t_PescanBase::escan.nlines    = nlines;
    t_PescanBase::escan.tolerance = tolerance;

    return result;
  }

  inline std::ostream& operator<<( std::ostream& _stream, const VirtualAtom& _va )
  {
    Ising_CE::Fourier( _va.structure.atoms.begin(),  _va.structure.atoms.end(),
                       _va.structure.k_vecs.begin(), _va.structure.k_vecs.end() );
    return _stream << _va.structure;
  }

} // namespace Pescan

#endif


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
   *               }_{\Delta \mathcal{H}_{\mathrm{stress}}}
   *               +
   *               \underbrace{
   *                 \sum_{\mathbf{R}}\frac{\partial \mathcal{H} }
   *                                       {\partial\epsilon_\mathcal{R}}
   *                                  \frac{\partial\epsilon_\mathcal{R}}
   *                                       {\partial S}
   *               }_{\Delta \mathcal{H}_{\mathrm{strain}}}
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
   *                   .
   *             </DD>
   *             <DT> \f$\Delta \mathcal{H}_{\mathrm{stress}}\f$
   *             <DD> 
   *                   - \b Keep the wavefunctions, unit-cell, and microscopic
   *                     strain from \f$sigma\f$
   *                   - \b Compute the strain relaxed Vff positions for
   *                     \f$\sigma'\f$ with constant unit-cell.
   *                   - \b Construct a perturbed potential consisting of the
   *                   displaced atoms. 
   *                   - \b Apply the wave-functions \e as are onto the
   *                   perturbed potential.  Eg escan with 0 iterations.
   *                   .
   *             </DD>
   *             <DT> \f$\Delta \mathcal{H}_{\mathrm{strain}}\f$
   *             <DD> 
   *                   - \b Keep the wavefunctions, fractional coordinates, and
   *                   unit-cell from \f$\sigma\f$
   *                   - \b Compute the strain relaxed Vff positions for
   *                     \f$\sigma'\f$ with constant unit-cell.
   *                   - \b Construct a perturbed potential consisting of the
   *                   difference between the original \f$\sigma\f$ potential
   *                   and a potential with the new microscopic strains.
   *                   - \b Apply the wave-functions \e as are onto the
   *                   perturbed potential.  Eg escan with 0 iterations.
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
  class VirtualAtom : protected Interface, public function::VirtualAtom
  {
     protected:
       //! Type from which the VA functional is derived
       typedef VirtualAtom t_VABase;
       //! Type of the pescan interface class
       typedef Interface t_PescanBase;
       //! \brief Type of the Valence Force Field Functional class
       //! \details It is fitted into a class which adds positional pertubation
       //!          stuff
       typedef Vff::PescanPosGrad t_Vff;

       //! Amplitude of the numerical derivative
       const static types::t_real deriv_amplitude = 0.01;

       //! The type of atom container
       typedef t_VABase :: t_Atoms t_Atoms;
       //! The type of atoms 
       typedef t_VABase :: t_Atom  t_Atom;

     protected:
       //! the vff object to get the variation of the positions
       t_Vff vff;


     public:
       //! Constructor and Initializer
       VirtualAtom   ( Ising_CE::Structure &_str )
                   : t_PescanBase(), t_VABase( _str ), 
                     vff( structure ) {}
       //! Copy Constructor
       VirtualAtom   ( const VirtualAtom &_c )
                   : t_PescanBase( _c ), t_VABase( _str ), 
                     vff(_c.structure) {}
        

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
       //! \brief partial derivative of the band gap with repect to occupation,
       //!        with constant positions.
       t_Type potential_grad( types :: t_int _pos );
  }


} // namespace Pescan

#endif


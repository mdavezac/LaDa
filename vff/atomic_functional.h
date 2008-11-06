//
//  Version: $Id$
//
#ifndef _VFF_ATOMIC_FUNCTIONAL_H_
#define _VFF_ATOMIC_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <algorithm>
#include <boost/filesystem/path.hpp>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

//! \brief Reimplements the Valence Force Field %Functional in c++
//! \details Vff, or Valence Force Field functional, is an empirical functional which
//! attempts to model the strain energy of a material from at most three-body
//! interactions. The there body interactions which are considered are
//! bond-stretching, change in bond-angles, and a combination of these two.
//! 
//! The implementation relies on a body-centered paradigm. In other words,
//! four classes have been created:
//!   - Vff::Functional is wrapper class and interface to Vff
//!   - Vff::Atomic_Center represent a single atom and lists its first neighbor relationships
//!   - Vff::Atomic_Center::const_iterator allows coders to travel
//!   through a collection of Vff::Atomic_Centera along first neighbor
//!   relationships.
//!   - Vff::Atomic_Functional computes the strain energy of one single atom, eg
//!   all three body terms in which a particular atom takes part.
//!   .
namespace Vff
{
  //! \cond
  class Atomic_Center;
  //! \endcond

  //! \brief An atom centered functional
  //! \details Atomic_Functional will compute the energy and strain resulting from
  //! bond-stretching, angle warping, and bond-angle-warping, eg from all
  //! first-neighbor two and three body interactions, on an Vff::Atomic_Center atom.
  class Atomic_Functional 
  {
    //! The type of the atom  
    typedef Crystal::Structure::t_Atom  t_Atom;
    const static types::t_real twos3;  //!<\f$2*\sqrt(3)\f$
    const static types::t_real one16;  //!<\f$\frac{1}{16}\f$
    const static types::t_real s3o160; //!<\f$\frac{\sqrt(3)}{8}\f$
    const static types::t_real one640; //!<\f$\frac{1}{640}\f$
    const static types::t_real three8; //!<\f$\frac{3}{8}\f$
    const static types::t_real s33o8;  //!<\f$\frac{3}{8}\sqrt(3)\f$
    const static types::t_real s33o16; //!<\f$\frac{3}{16}\sqrt(3)\f$
    const static types::t_real thre16; //!<\f$\frac{3}{16}\f$
    const static types::t_real thre32; //!<\f$\frac{3}{32}\f$
    const static types::t_real s33128; //!<\f$\frac{3}{128}\sqrt(3)\f$
    const static types::t_real s33256; //!<\f$\frac{3}{256}\sqrt(3)\f$
    const static types::t_real no1280; //!< some number, but which?
    const static types::t_real no2560; //!< some number, but which?

    protected:
      std::string specie;                  //!< atomic type as a string
      Crystal :: Structure *structure; //!< structure to which the Atomic_Functional belongs
      types::t_unsigned site;           //!< site number of Atomic_Center in Crystal::lattice
      types::t_unsigned type;           //!< atomic type of Atomic_Center
      std::vector< types::t_real > lengths;  //!< equilibrium bond lengths
      std::vector< types::t_real > alphas;   //!< bond stretching parameters
      std::vector< types::t_real > betas;    //!< angle deformation parameters
      std::vector< types::t_real > gammas;   //!< bond-angle parameters
      std::vector< types::t_real > sigmas;   //!< equilibrium tetrahedral symmetry
      
    public:
      //! Constructore.
      Atomic_Functional () : structure( NULL ) {};
      //! \brief Constructor and Initializer
      //! \param _str atomic type as a string
      //! \param _struct structure to which Atomic_Center belongs
      //! \param _site site number of Atomic_Center
      //! \param _type type number of Atomic_Center
      Atomic_Functional   ( std::string _str, Crystal::Structure &_struct, 
                            types::t_unsigned _site, 
                            types::t_unsigned _type )
                        : structure(&_struct), site(_site), type(_type) {}
      //! Copy Constructor
      Atomic_Functional   ( const Atomic_Functional &_a )
                        : specie(_a.specie), structure(_a.structure), site(_a.site), type(_a.type),
                          lengths(_a.lengths), alphas(_a.alphas),
                          betas(_a.betas), gammas(_a.gammas), sigmas(_a.sigmas) {}
      
      //! \brief Adds a bond type to bond list
      //! \param _typeB type of atom at end-point of bond
      //! \param _l equilibrium length
      //! \param _i bond stretching parameters
      void add_bond( const types::t_unsigned _typeB, const types::t_real _l,
                     const std::vector<types::t_real> &_i );
      //! \brief Adds a bond type to bond list
      //! \param _typeB type of atom at end-point of bond
      //! \param _l equilibrium length
      //! \param _i bond array of stretching parameters (5 types::t_real long)
      void add_bond( const types::t_unsigned _typeB, const types::t_real _l,
                     const types::t_real _i[5] );
      //! \brief Adds angle deformation and bond-angle parameters
      //! \param _typeA type of atom at end-point of one bond
      //! \param _typeC type of atom at end-point of other bond
      //! \param _gamma bond-angle deformation parameter
      //! \param _sigma equilibrium angle
      //! \param _i angle deformation parameters
      void add_angle( const types::t_unsigned _typeA,
                      const types::t_unsigned _typeC,
                      const types::t_real _gamma, const types::t_real _sigma, 
                      const std::vector<types::t_real> &_i );
      //! \brief Adds angle deformation and bond-angle parameters
      //! \param _typeA type of atom at end-point of one bond
      //! \param _typeC type of atom at end-point of other bond
      //! \param _gamma bond-angle deformation parameter
      //! \param _sigma equilibrium angle
      //! \param _i angle array of deformation parameters (5 types::t_real long)
      void add_angle( const types::t_unsigned _typeA,
                      const types::t_unsigned _typeC,
                      const types::t_real _gamma, const types::t_real _sigma, 
                      const types::t_real _i[5] );
      
      //! \brief Evaluate strain energy for Atomic_Center _center
      //! \details returns strain energy 
      //! \param _center center for which to evaluate energy
      //! \sa function::Base, function::Base::evaluate()
      types::t_real evaluate( const Atomic_Center &_center ) const;
      //! \brief Evaluate strain energy and gradients for Atomic_Center _center
      //! \details returns strain energy, and computes stress
      //! \param _center center for which to evaluate energy
      //! \param _strain to which Atomic_Functional::structure is submitted
      //! \param _stress on _center resulting from _strain
      //! \param _K0 displacement resulting from _strain and w.r.t original unit-cell
      //! \sa function::Base, function::Base::evaluate_with_gradient()
      types::t_real evaluate_with_gradient( Atomic_Center &_center,
                                            const atat::rMatrix3d &_strain,
                                            atat::rMatrix3d &_stress,
                                            const atat::rMatrix3d &_K0 ) const;
      //! \brief computes the trace of the microscopic strain on an atomic center
      //!  structure0 and the atomic centers are expected to be related 
      //! \details To be used for pescan
      types::t_real MicroStrain( const Atomic_Center &_center, 
                                 const Crystal::Structure &_str0 ) const;

      //! prints out all parameters
      void print_out( std::ostream &stream ) const;

      //! \brief Loads an atomic center from XML.
      //! \param _node is \<Functional type=\"vff\"\>.
      //! \param _site_index is the index of the lattice-site for this atomic functional.
      //! \param _site_type is the type of the lattice-site for this atomic functional.
      //! \param _othersite is the lattice-site of the first neighbours.
      bool load( const TiXmlElement& _node,
                 const types::t_unsigned &_site_index,
                 const types::t_unsigned &_type_index,
                 const Crystal::Lattice::t_Site &_othersite 
                 Crystal :: Structure &_structure );
  }; 
} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_

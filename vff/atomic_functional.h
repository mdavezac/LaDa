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

namespace LaDa
{
  namespace Vff
  {
    //! \cond
    class AtomicCenter;
    //! \endcond

    //! \brief An atom centered functional
    //! \details AtomicFunctional will compute the energy and strain resulting from
    //! bond-stretching, angle warping, and bond-angle-warping, eg from all
    //! first-neighbor two and three body interactions, on an Vff::AtomicCenter atom.
    class AtomicFunctional 
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
        Crystal :: Structure *structure; //!< structure to which the AtomicFunctional belongs
        types::t_unsigned site;           //!< site number of AtomicCenter in Crystal::lattice
        types::t_unsigned type;           //!< atomic type of AtomicCenter
        std::vector< types::t_real > lengths;  //!< equilibrium bond lengths
        std::vector< types::t_real > alphas;   //!< bond stretching parameters
        std::vector< types::t_real > betas;    //!< angle deformation parameters
        std::vector< types::t_real > gammas;   //!< bond-angle parameters
        std::vector< types::t_real > sigmas;   //!< equilibrium tetrahedral symmetry
        
      public:
        //! Constructore.
        AtomicFunctional () : structure( NULL ) {};
        //! \brief Constructor and Initializer
        //! \param _str atomic type as a string
        //! \param _struct structure to which AtomicCenter belongs
        //! \param _site site number of AtomicCenter
        //! \param _type type number of AtomicCenter
        AtomicFunctional   ( std::string _str, Crystal::Structure &_struct, 
                              types::t_unsigned _site, 
                              types::t_unsigned _type )
                          : structure(&_struct), site(_site), type(_type) {}
        //! Copy Constructor
        AtomicFunctional   ( const AtomicFunctional &_a )
                          : specie(_a.specie), structure(_a.structure), site(_a.site), type(_a.type),
                            lengths(_a.lengths), alphas(_a.alphas),
                            betas(_a.betas), gammas(_a.gammas), sigmas(_a.sigmas) {}
        
        //! \brief Evaluate strain energy for AtomicCenter _center
        //! \details returns strain energy 
        //! \param _center center for which to evaluate energy
        //! \sa function::Base, function::Base::evaluate()
        types::t_real evaluate( const AtomicCenter &_center ) const;
        //! \brief Evaluate strain energy and gradients for AtomicCenter _center
        //! \details returns strain energy, and computes stress
        //! \param _center center for which to evaluate energy
        //! \param _strain to which AtomicFunctional::structure is submitted
        //! \param _stress on _center resulting from _strain
        //! \param _K0 displacement resulting from _strain and w.r.t original unit-cell
        //! \sa function::Base, function::Base::evaluate_with_gradient()
        types::t_real evaluate_with_gradient( const AtomicCenter &_center,
                                              const atat::rMatrix3d &_strain,
                                              atat::rMatrix3d &_stress,
                                              const atat::rMatrix3d &_K0 ) const;
        //! \brief computes the trace of the microscopic strain on an atomic center
        //!  structure0 and the atomic centers are expected to be related 
        //! \details To be used for pescan
        types::t_real MicroStrain( const AtomicCenter &_center, 
                                   const Crystal::Structure &_str0 ) const;

        //! prints out all parameters
        void print_out( std::ostream &stream ) const;

        //! \brief Loads an atomic center from XML.
        //! \param _node is \<Functional type=\"vff\"\>.
        //! \param _site_index is the index of the lattice-site for this atomic functional.
        //! \param _type_index is the type of the lattice-site for this atomic functional.
        //! \param _othersite is the lattice-site of the first neighbours.
        //! \param _structure Needed for its scale only. The AtomicFunctional
        //!                   instance will keep a pointer to this object.
        bool load( const TiXmlElement& _node,
                   const types::t_unsigned &_site_index,
                   const types::t_unsigned &_type_index,
                   const Crystal::Lattice::t_Site &_othersite,
                   Crystal :: Structure &_structure );
    }; 
  } // namespace vff 
} // namespace LaDa
#endif // _VFF_FUNCTIONAL_H_

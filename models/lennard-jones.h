//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _MODELS_LENNARDJONES_FUNCTIONAL_H_
#define _MODELS_LENNARDJONES_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <opt/debug.h>


namespace LaDa
{
  namespace Models
  {
    //! \brief A LennardJones model functional.
    //! \details Ad-hoc implementation for Clj functional.
    class LennardJones
    {
      public:
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;
        //! Constructor and Initializer
        LennardJones() {}
        //! Copy Constructor
        LennardJones   ( const LennardJones &_c )
                     : species_(_c.species_), bond_strength_(_c.bond_strength_),
                       mesh_(_c.mesh_), rcut_(_c.rcut_) {}
        //! \brief Destructor
        LennardJones() {}

        //! \brief computes energy and stress and forces.
        //! \param[in] _in input structure, with reduced coordinates.
        //! \param[inout] _out stress and forces. Both are \e added to existing stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

      protected:
        //! A structure holding the specie types.
        struct Bond;
        //! Type of the key used in the map.
        typedef std::string Key;
        //! Type of the container of atomic species.
        typedef std::map< Key, Bond > t_Bonds;

        //! Contains all atomic species.
        t_Bonds species_;
        //! Cutoff mesh.
        atat::rVector3d mesh_;
        //! Real space cutoff.
        types::t_real rcut_;
    };

    struct LennardJones :: Bond
    {
      //! Type of the atomic mass.
      typedef std::string t_Type;
      //! Atomic Mass.
      t_Type type;
      //! Hard sphere parameter.
      types::t_real hard_sphere;
      //! van der Walls parameter.
      types::t_real vand_der_walls;
      //! Constructor.
      Bond   ( const t_Type &_t, const types::t_real &_hs, const types::t_real &_vdw )
             : type( _t ), hard_sphere( _r ), van_der_walls( _vdw ) {}
      //! Copy Constructor.
      Bond   ( const Bond& _t )
             : type( _c.type ), hard_sphere( _c.hard_sphere ), van_der_walls( _c.van_der_walls ) {}
    };

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_

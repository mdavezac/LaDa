//
//  Version: $Id: vff.h 895 2008-12-22 02:04:18Z davezac $
//
#ifndef _MODELS_LENNARDJONES_FUNCTIONAL_H_
#define _MODELS_LENNARDJONES_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <boost/filesystem/path.hpp>
#include <boost/tuple/tuple.hpp>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>


namespace LaDa
{
  namespace Models
  {
    //! \brief A LennardJones model functional.
    //! \details Ad-hoc implementation for Clj functional.
    class LennardJones
    {
      protected:
        //! Cutoff mesh type.
        typedef boost::tuple<types::t_int, types::t_int, types::t_int > t_MeshTuple;

      public:
        //! A structure holding the specie types.
        struct Bond;
        //! Type of the key used in the map.
        typedef std::string Key;
        //! Type of the container of atomic species.
        typedef std::map< Key, Bond > t_Bonds;
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;


        //! Contains all atomic species.
        t_Bonds bonds;

        //! Constructor and Initializer
        LennardJones() {}
        //! Copy Constructor
        LennardJones   ( const LennardJones &_c )
                     : bonds(_c.bonds), mesh_(_c.mesh_), rcut_(_c.rcut_) {}
        //! Destructor.
        ~LennardJones() {}

        //! \brief computes energy and stress and forces.
        //! \param[in] _in input structure, with reduced coordinates.
        //! \param[inout] _out stress and forces. Both are \e added to existing stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

        //! Loads parameters from XML.
        bool Load( const TiXmlElement& _node );

        //! Gets mesh for cutoff.
        t_MeshTuple get_mesh() const { return mesh_; }
        //! Sets mesh for cutoff.
        void set_mesh( const t_MeshTuple &_t ) { mesh_ = _t; }
        //! Gets cutoff.
        types::t_real get_rcutoff() const { return rcut_; }
        //! Sets cutoff.
        void set_rcutoff(types::t_real _r) { rcut_ = _r; }


      protected:

        //! Returns bond name from bond endpoints.
        static Key bondname( const Key &_a, const Key &_b )
          { return _a > _b ? _a + " " +  _b: _b + " " + _a; }
        //! Extracts atom A from bond name.
        Key extract_atomA( const Key &_key ) const
          { return _key.substr( 0, _key.find( " " ) ); }
        //! Extracts atom B from bond name.
        Key extract_atomB( const Key &_key ) const
          { return _key.substr( _key.find( " " ) ); }

        //! Makes sure that bonds are coherently declared.
        void check_coherency() const;

        //! Cutoff mesh.
        t_MeshTuple mesh_;
        //! Real space cutoff.
        types::t_real rcut_;
    };

    struct LennardJones :: Bond
    {
      //! Hard sphere parameter.
      types::t_real hard_sphere;
      //! van der Walls parameter.
      types::t_real van_der_walls;
      //! Constructor.
      Bond() {};
      //! Constructor.
      Bond   ( const types::t_real &_hs, const types::t_real &_vdw )
             : hard_sphere( _hs ), van_der_walls( _vdw ) {}
      //! Copy Constructor.
      Bond   ( const Bond& _c )
             : hard_sphere( _c.hard_sphere ), van_der_walls( _c.van_der_walls ) {}
      //! Loads parameters from XML.
      bool Load( const TiXmlElement& _node, std::string &_type );
    };

  } // namespace CLJ.
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_

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
#include "python/clj.hpp"


namespace LaDa
{
  namespace Models
  {
    //! \cond
    class Clj;
    void read_fortran_input( Clj&, const boost::filesystem::path& );
    //! \endcond

    //! \brief A LennardJones model functional.
    //! \details Ad-hoc implementation for Clj functional.
    class LennardJones
    {
      FRIEND_EXPOSE_CLJ
      friend  void read_fortran_input( Clj&, const boost::filesystem::path&);
      public:
        //! Argument type.
        typedef Crystal :: TStructure< std::string > t_Arg;
        //! Return type.
        typedef types::t_real t_Return;
        //! Constructor and Initializer
        LennardJones() {}
        //! Copy Constructor
        LennardJones   ( const LennardJones &_c )
                     : bonds_(_c.bonds_), mesh_(_c.mesh_), rcut_(_c.rcut_) {}
        //! Destructor.
        ~LennardJones() {}

        //! \brief computes energy and stress and forces.
        //! \param[in] _in input structure, with reduced coordinates.
        //! \param[inout] _out stress and forces. Both are \e added to existing stress and forces.
        t_Return energy(const t_Arg& _in, t_Arg& _out) const;

        //! Loads parameters from XML.
        bool Load( const TiXmlElement& _node );

      protected:
        //! Type of the key used in the map.
        typedef std::string Key;

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

        //! A structure holding the specie types.
        struct Bond;
        //! Type of the container of atomic species.
        typedef std::map< Key, Bond > t_Bonds;

        //! Contains all atomic species.
        t_Bonds bonds_;
        //! Cutoff mesh type.
        typedef boost::tuple<types::t_int, types::t_int, types::t_int > t_MeshTuple;
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

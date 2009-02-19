//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include "ewald.h"

extern "C" void FC_FUNC( ewaldf, EWALDF )
                (
                  const double *const, // verbosity
                  const double *const, // Energy
                  const double *const, // forces (reduced)
                  const double *const, // forces (cartesian)
                  const double *const, // stress
                  const double *const, // number of atoms
                  const double *const, // reduced atomic coordinates.
                  const double *const, // atomic charges
                  const double *const, // cell vectors
                  const double *const  // dimension of arrays.
                );
namespace LaDa
{
  namespace Models
  {
    Ewald :: t_Return Ewald :: energy( const t_Arg& _in, t_Arg &_out ) const
    {
      __DOASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )

      const size_t natoms( _in.atoms.size() );
      double charges[ natoms ];
      double positions[ natoms * 3 ], forces[ natoms * 3 ], cforces[ natoms * 3 ];
      double cell[ 9 ], stress[ 6 ];
      const double verbosity(0);
      double energy(0);
      const double n( natoms );

      typedef t_Arg :: t_Atoms :: const_iterator t_cit;
      t_cit i_atom = _in.atoms.begin();
      for( size_t i(0); i < natoms; ++i, ++i_atom )
      {
        __ASSERT( charges_.find( i_atom->type ) == charges_.end(),
                  "Atomic charge does not exist.\n" )
        charges[i] = charges_.find(i_atom->type)->second;
        positions[ i*3 ]     = i_atom->pos[0];
        positions[ i*3 + 1 ] = i_atom->pos[1];
        positions[ i*3 + 2 ] = i_atom->pos[2];
      } // loop over atoms

      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          cell[ i + j*3 ] = _in.cell(i,j);

      FC_FUNC( ewaldf, EWALDF )
      (
        &verbosity,   // verbosity
        &energy,      // Energy
        forces,       // forces (reduced)
        cforces,      // forces (cartesian)
        stress,       // stress
        &n,           // number of atoms
        positions,    // reduced atomic coordinates.
        charges,      // atomic charges
        cell,         // cell vectors
        &n            // dimension of arrays.
      );
      // Copy (reduced) forces.
      t_Arg :: t_Atoms :: iterator i_force = _out.atoms.begin();
      for( size_t i(0); i < natoms; ++i, ++i_force )
      {
        i_force->pos[0] += forces[ i*3 ];
        i_force->pos[1] += forces[ i*3 + 1 ];
        i_force->pos[2] += forces[ i*3 + 2 ];
      } // loop over atoms
      // copy stress.
      _out.cell(0,0) += stress[0];
      _out.cell(1,1) += stress[1];
      _out.cell(2,2) += stress[2];
      _out.cell(0,1) += stress[3];
      _out.cell(1,0) += stress[3];
      _out.cell(1,2) += stress[4];
      _out.cell(2,1) += stress[4];
      _out.cell(0,2) += stress[5];
      _out.cell(2,0) += stress[5];
      
      return energy;
    }

    bool Ewald :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Ewald" );
      if( not parent ) return false;
      const TiXmlElement* child = parent->FirstChildElement( "Atom" );
      for(; child; child = child->NextSiblingElement("Atom") )
      {
        if( child->Attribute("type") and child->Attribute("Charge") ) continue;
        const std::string type = Print::StripEdges( child->Attribute("type") );
        const types::t_real charge
          = boost::lexical_cast<types::t_real>( child->Attribute("charge") );
        __DOASSERT( charges_.find( type ) != charges_.end(),
                    "Duplicate entry in Ewald functional for " + type + "\n" )
        charges_[type] = charge;
      }

      return true;
    }

  } // namespace Models
} // namespace LaDa

#include "LaDaConfig.h"
#include "FCMangle.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <physics/physics.h>

#include "ewald.h"

extern "C" void FC_GLOBAL( ewaldf, EWALDF )
                (
                  const int *const, // verbosity
                  const double *const, // Energy
                  const double *const, // forces (reduced)
                  const double *const, // forces (cartesian)
                  const double *const, // stress
                  const int *const, // number of atoms
                  const double *const, // reduced atomic coordinates.
                  const double *const, // atomic charges
                  const double *const, // real space cutoff
                  const double *const, // cell vectors
                  const int *const  // dimension of arrays.
                );
namespace LaDa
{
  namespace Models
  {
    Ewald :: t_Return Ewald :: operator()( const t_Arg& _in, t_Arg &_out ) const
    {
      LADA_DO_NASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )

      const size_t natoms( _in.atoms.size() );
      double Qs[ natoms ];
      double positions[ natoms * 3 ], forces[ natoms * 3 ], cforces[ natoms * 3 ];
      double cell[ 9 ], stress[ 6 ];
      const int verbosity(0);
      double energy(0);
      const int n( natoms );
      const double rcut( cutoff_ / Physics::a0("A") );
      math::rMatrix3d const inv(!_in.cell);

      typedef t_Arg :: t_Atoms :: const_iterator t_cit;
      t_cit i_atom = _in.atoms.begin();
      for( size_t i(0); i < natoms; ++i, ++i_atom )
      {
        LADA_NASSERT( charges.find( i_atom->type ) == charges.end(),
                  "Atomic charge does not exist.\n" )
        Qs[i] = charges.find(i_atom->type)->second;
        math::rVector3d const frac(inv * i_atom->pos);
        positions[ i*3 ]     = frac[0];
        positions[ i*3 + 1 ] = frac[1];
        positions[ i*3 + 2 ] = frac[2];
      } // loop over atoms

      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          cell[ i + j*3 ] = _in.cell(i,j) * _in.scale / Physics::a0("A");

      FC_GLOBAL( ewaldf, EWALDF )
      (
        &verbosity,   // verbosity
        &energy,      // Energy
        forces,       // forces (reduced)
        cforces,      // forces (cartesian)
        stress,       // stress
        &n,           // number of atoms
        positions,    // reduced atomic coordinates.
        Qs,           // atomic charges
        &rcut,        // Real-space cutoff.
        cell,         // cell vectors
        &n            // dimension of arrays.
      );
      // Copy (reduced) forces.
      t_Arg :: t_Atoms :: iterator i_force = _out.atoms.begin();
      for( size_t i(0); i < natoms; ++i, ++i_force )
      {
        i_force->pos[0] += cforces[ i*3 ]     * Physics::Rydberg("eV") / Physics::a0("A");
        i_force->pos[1] += cforces[ i*3 + 1 ] * Physics::Rydberg("eV") / Physics::a0("A");
        i_force->pos[2] += cforces[ i*3 + 2 ] * Physics::Rydberg("eV") / Physics::a0("A");
      } // loop over atoms
      // copy stress.
      _out.scale = 1e0;
      _out.cell(0,0) += stress[0] * Physics::Rydberg("eV");
      _out.cell(1,1) += stress[1] * Physics::Rydberg("eV");
      _out.cell(2,2) += stress[2] * Physics::Rydberg("eV");
      _out.cell(0,1) += stress[3] * Physics::Rydberg("eV");
      _out.cell(1,0) += stress[3] * Physics::Rydberg("eV");
      _out.cell(1,2) += stress[4] * Physics::Rydberg("eV");
      _out.cell(2,1) += stress[4] * Physics::Rydberg("eV");
      _out.cell(0,2) += stress[5] * Physics::Rydberg("eV");
      _out.cell(2,0) += stress[5] * Physics::Rydberg("eV");
      
      return energy * Physics::Rydberg("eV");
    }

    bool Ewald :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Ewald" );
      if( not parent ) return false;
      const TiXmlElement* child = parent->FirstChildElement( "Atom" );
      for(; child; child = child->NextSiblingElement("Atom") )
      {
        if( child->Attribute("type") and child->Attribute("Charge") ) continue;
        const std::string type = boost::algorithm::trim_copy(std::string(child->Attribute("type")));
        const types::t_real charge
          = boost::lexical_cast<types::t_real>( child->Attribute("charge") );
        LADA_DO_NASSERT( charges.find( type ) != charges.end(),
                    "Duplicate entry in Ewald functional for " + type + "\n" )
        charges[type] = charge;
      }

      return true;
    }

  } // namespace Models
} // namespace LaDa

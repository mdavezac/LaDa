#include "LaDaConfig.h"

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <opt/tuple_io.h>

#include "lennard-jones.h"

namespace LaDa
{
  namespace Models
  {
    namespace details
    {
      types::t_real range01( const types::t_real _a )
      {
        const types::t_real b( _a - std::floor( _a ) );
        return _a > 1e0 ? _a - 1e0: _a < 0e0 ? _a + 1e0: _a; 
      }
    }
    LennardJones :: t_Return LennardJones :: operator()( const t_Arg& _in, t_Arg &_out ) const
    {
      namespace bt = boost::tuples;
      LADA_DO_NASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )
      types::t_real energy(0), ecut(0);
      const types::t_real scale_squared( _in.scale * _in.scale );

      typedef t_Arg ::  t_Atoms :: const_iterator t_cit;
      const t_cit i_atom_begin = _in.atoms.begin();
      const t_cit i_atom_end = _in.atoms.end();
      const t_Arg :: t_Atoms :: iterator i_force_begin = _out.atoms.begin();
      t_Arg :: t_Atoms :: iterator i_force1( i_force_begin );
      for( t_cit i_atom1( i_atom_begin ); i_atom1 != i_atom_end; ++i_atom1, ++i_force1 )
      {
        t_Arg :: t_Atoms :: iterator i_force2(i_force_begin);
        for( t_cit i_atom2( i_atom_begin ); i_atom2 != i_atom_end; ++i_atom2, ++i_force2 )
        {
          const std::string bondtype( bondname( i_atom2->type, i_atom1->type ) );

          LADA_DO_NASSERT( bonds.end() == bonds.find( bondtype ),
                      "Bond " + i_atom1->type + "-" + i_atom2->type + " does not exist.\n" )
          const Bond &bond( bonds.find( bondtype )->second );
          const types::t_real rcut_squared( rcut_ * rcut_ );

          const math::rVector3d dfractional
                                ( 
                                  details::range01( i_atom1->pos[0] - i_atom2->pos[0] ),
                                  details::range01( i_atom1->pos[1] - i_atom2->pos[1] ),
                                  details::range01( i_atom1->pos[2] - i_atom2->pos[2] )
                                );
          for( types::t_int i(-bt::get<0>(mesh_)); i < bt::get<0>(mesh_); ++i )
            for( types::t_int j(-bt::get<1>(mesh_)); j < bt::get<1>(mesh_); ++j )
              for( types::t_int k(-bt::get<2>(mesh_)); k < bt::get<2>(mesh_); ++k )
              {
                if ( i_atom1 == i_atom2 and i == 0 and j == 0 and k == 0 ) continue;
                // computes distance.
                const math::rVector3d fdistance( dfractional + math::rVector3d(i,j,k) );
                const math::rVector3d distance( _in.cell * fdistance * _in.scale );
                const types::t_real normd( distance.squaredNorm() );
                if( normd > rcut_squared ) continue;

                // result -= 4.0 * scale * lj_pot( sigma_squared / rcut_squared )
                { // compute cutoff correction energy
                  const types::t_real squared( 1e0 / rcut_squared );
                  const types::t_real sixth( squared * squared * squared );
                  const types::t_real twelfth( sixth * sixth );
                  ecut -=  bond.hard_sphere * twelfth - bond.van_der_walls * sixth;
                }

                { // van_der_walls energy, force, stress.
                  const types::t_real squared( 1e0 / normd );
                  const types::t_real sixth( squared * squared * squared );
                  const types::t_real twelfth( sixth * sixth );
                  energy += bond.hard_sphere * twelfth - bond.van_der_walls * sixth;

                  const types::t_real ffactor
                  ( 
                      6e0 * squared
                    * ( 2.e0 * bond.hard_sphere * twelfth - bond.van_der_walls * sixth )
                  );
                  const math::rVector3d force( ffactor * distance );
                  i_force1->pos += force;
                  i_force2->pos -= force;

                  for( size_t i(0); i < 3; ++i )
                    for( size_t j(0); j < 3; ++j )
                      _out.cell(i,j) += force[i] * distance[j];
                }

              } // loop over periodic images.
        } // loop over atom 2
      } // loop over atom 1
      return energy + ecut;
    }
    
    bool LennardJones :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "LennardJones" );
      LADA_DO_NASSERT( not parent->Attribute( "mesh" ), 
                  "LennardJones functional requires a mesh attribute.\n" )
      LADA_DO_NASSERT( not parent->Attribute( "rcut" ), 
                  "LennardJones functional requires a rcut attribute.\n" )

      // reads bonds
      const TiXmlElement* child = opt::find_node( _node, "Bond");
      bonds.clear();
      for(; child; child = child->NextSiblingElement( "Bond" ) )
      {
        Bond bond;
        std::string type("");
        if( not bond.Load( *child, type ) ) return false;
        LADA_DO_NASSERT( bonds.end() != bonds.find( type ),
                    "Duplicate bond type.\n" )
        bonds[type] = bond;
      }
       
      // reads attributes.
      LADA_DO_NASSERT( not opt::tuples::read( parent->Attribute("mesh"), mesh_ ),
                  "Could not parse mesh attribute.\n" )
      rcut_ = boost::lexical_cast< types::t_real >( parent->Attribute("rcut") );
      return true;
    }

    bool LennardJones :: Bond :: Load( const TiXmlElement& _node, std::string &_type )
    {
      LADA_NASSERT( _node.Value() != "Bond", "Incorrect XML node.\n" )
      LADA_DO_NASSERT( _node.Attribute( "A" ), "Bond requires an A attribute.\n" )
      LADA_DO_NASSERT( _node.Attribute( "B" ), "Bond requires a B attribute.\n" )
      LADA_DO_NASSERT( _node.Attribute( "hardsphere" ), "Bond requires a hardsphere attribute.\n" )
      LADA_DO_NASSERT( _node.Attribute( "vanderwalls" ), "Bond requires a vanderwalls attribute.\n" )
      const std::string A = boost::algorithm::trim_copy(std::string(_node.Attribute("A")));
      const std::string B = boost::algorithm::trim_copy(std::string(_node.Attribute("B")));
      _type = bondname(A, B);
      hard_sphere = boost::lexical_cast< types::t_real >( _node.Attribute("hard_sphere") );
      van_der_walls = boost::lexical_cast< types::t_real >( _node.Attribute("van_der_walls") );
    }


    void LennardJones :: check_coherency() const
    {
      std::vector< std::string > lj;
      t_Bonds :: const_iterator i_bond = bonds.begin();
      t_Bonds :: const_iterator i_bond_end = bonds.end();
      for(; i_bond != i_bond_end; ++i_bond) 
      {
        const std::string A( LennardJones :: extract_atomA( i_bond->first ) );
        const std::string B( LennardJones :: extract_atomB( i_bond->first ) );
        if( lj.end() == std::find( lj.begin(), lj.end(), A ) ) lj.push_back( A );
        if( lj.end() == std::find( lj.begin(), lj.end(), B ) ) lj.push_back( B );
      }
      bool result = true;
      foreach( const Key &A, lj )
        foreach( const Key &B, lj )
        {
          if( bonds.find( bondname(A, B) ) != bonds.end() ) continue;
          std::cerr << "Bond " + bondname( A, B ) + " does not exist.\n";
          result = false;
        }

      LADA_DO_NASSERT( not result, "" );
    }

  } // namespace CLJ
} // namespace LaDa

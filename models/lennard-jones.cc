//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include <print/manip.h>
#include <opt/tuple_io.h>

#include "lennard-jones.h"

namespace LaDa
{
  namespace Models
  {
    LennardJones :: t_Return LennardJones :: energy( const t_Arg& _in, t_Arg &_out ) const
    {
      namespace bt = boost::tuples;
      __DOASSERT( _in.atoms.size() != _out.atoms.size(), "Incoherent structure size.\n" )
      types::t_real energy(0);

      typedef t_Arg ::  t_Atoms :: const_iterator t_cit;
      t_cit i_atom_begin = _in.atoms.begin();
      t_cit i_atom_end = _in.atoms.end();
      t_Arg :: t_Atoms :: iterator i_force1 = _out.atoms.begin();
      for( t_cit i_atom1( i_atom_begin ); i_atom1 != i_atom_end; ++i_atom1, ++i_force1 )
      {
        t_Arg :: t_Atoms :: iterator i_force2 = i_force1 + 1;
        for( t_cit i_atom2( i_atom1 + 1 ); i_atom2 != i_atom_end; ++i_atom2, ++i_force2 )
        {
          const std::string bondtype( i_atom2->type > i_atom1->type ?
                                      i_atom2->type + i_atom1->type:
                                      i_atom1->type + i_atom2->type );

          __DOASSERT( bonds_.end() != bonds_.find( bondtype ),
                      "Bond " + i_atom1->type + "-" + i_atom2->type + " does not exist.\n" )
          const Bond &bond( bonds_.find( bondtype )->second );
          const types::t_real rcut_squared( rcut_ * rcut_ );

          const atat::rVector3d dfractional
                                ( 
                                  std::floor( i_atom1->pos[0] - i_atom2->pos[0] ),
                                  std::floor( i_atom1->pos[1] - i_atom2->pos[1] ),
                                  std::floor( i_atom1->pos[2] - i_atom2->pos[2] )
                                );
          for( size_t i(-bt::get<0>(mesh_)); i < bt::get<0>(mesh_); ++i )
            for( size_t j(-bt::get<1>(mesh_)); j < bt::get<1>(mesh_); ++j )
              for( size_t k(-bt::get<2>(mesh_)); k < bt::get<2>(mesh_); ++k )
              {
                // computes distance.
                const atat::rVector3d distance
                (
                  _in.cell * ( dfractional + atat::rVector3d(i,j,k) ) 
                );
                const types::t_real normd( atat::norm2(distance) );
                if( normd > rcut_squared ) continue;

                // energy -= 4.0 * scale * lj_pot( sigma_squared / rcut_squared )
                { // compute cutoff correction energy
                  const types::t_real squared( 1e0 / rcut_squared );
                  const types::t_real sixth( squared * squared * squared );
                  const types::t_real twelveth( sixth * sixth );
                  energy -=  bond.hard_sphere * twelveth - bond.van_der_walls * sixth;
                }

                { // van_der_walls energy, force, stress.
                  const types::t_real squared( 1e0 / normd );
                  const types::t_real sixth( squared * squared * squared );
                  const types::t_real twelveth( sixth * sixth );
                  energy -= bond.hard_sphere * twelveth - bond.van_der_walls * sixth;

                  const types::t_real ffactor
                  ( 
                    squared * ( 2.e0 * bond.hard_sphere * twelveth - bond.van_der_walls * sixth )
                  );
                  const atat::rVector3d force( ffactor * distance );
                  i_force1->pos -= force;
                  i_force2->pos += force;

                  for( size_t i(0); i < 3; ++i )
                    for( size_t j(0); j < 3; ++j )
                      _out.cell(i,j) += -force[i] * distance[j];
                }

              } // loop over periodic images.
        } // loop over atom 2
      } // loop over atom 1
    }
    
    bool LennardJones :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "LennardJones" );
      __DOASSERT( not parent->Attribute( "mesh" ), 
                  "LennardJones functional requires a mesh attribute.\n" )
      __DOASSERT( not parent->Attribute( "rcut" ), 
                  "LennardJones functional requires a rcut attribute.\n" )

      // reads bonds
      const TiXmlElement* child = opt::find_node( _node, "Bond");
      bonds_.clear();
      for(; child; child = child->NextSiblingElement( "Bond" ) )
      {
        Bond bond;
        if( not bond.Load( *child ) ) return false;
        __DOASSERT( bonds_.end() != bonds_.find( bond.type ),
                    "Duplicate bond type.\n" )
        bonds_[bond.type] = bond;
      }
       
      // reads attributes.
      __DOASSERT( not opt::tuples::read( parent->Attribute("mesh"), mesh_ ),
                  "Could not parse mesh attribute.\n" )
      rcut_ = boost::lexical_cast< types::t_real >( parent->Attribute("rcut") );
      return true;
    }

    bool LennardJones :: Bond :: Load( const TiXmlElement& _node )
    {
      __ASSERT( _node.Value() != "Bond", "Incorrect XML node.\n" )
      __DOASSERT( _node.Attribute( "A" ), "Bond requires an A attribute.\n" )
      __DOASSERT( _node.Attribute( "B" ), "Bond requires a B attribute.\n" )
      __DOASSERT( _node.Attribute( "hardsphere" ), "Bond requires a hardsphere attribute.\n" )
      __DOASSERT( _node.Attribute( "vanderwalls" ), "Bond requires a vanderwalls attribute.\n" )
      const std::string A = Print :: StripEdges( _node.Attribute("A") );
      const std::string B = Print :: StripEdges( _node.Attribute("B") );
      type = A > B ? A + B: B + A;
      hard_sphere = boost::lexical_cast< types::t_real >( _node.Attribute("hard_sphere") );
      van_der_walls = boost::lexical_cast< types::t_real >( _node.Attribute("van_der_walls") );
    }

  } // namespace CLJ
} // namespace LaDa

//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/structure.h>
#include <crystal/neighbors.h>

#include "bases.h"
#include "representation.h"

namespace LaDa
{
  namespace atomic_potential
  {
    // Strict weak ordering functor from knowledge of basis.
    struct basis_sort
    {
      Basis const& basis;
      basis_sort( Basis const &_b ) : basis(_b) {}
      basis_sort( basis_sort const &_b ) : basis(_b.basis) {}
      bool operator()( Crystal::Neighbor const * const _a,
                       Crystal::Neighbor const * const _b ) const
      {
        if( Fuzzy::neq( _a->distance, _b->distance ) ) return _a->distance < _b->distance;
        const types::t_real ax = _a->pos * basis.x;
        const types::t_real bx = _b->pos * basis.x;
        if( Fuzzy::neq( ax, bx ) ) return ax > bx;
        const types::t_real ay = _a->pos * basis.y;
        const types::t_real by = _b->pos * basis.y;
        if( Fuzzy::neq( ay, by ) ) return ay > by;
        const types::t_real az = _a->pos * basis.z;
        const types::t_real bz = _b->pos * basis.z;
        if( Fuzzy::neq( az, bz ) ) return az > bz;
        return false;
      }
    }; 


    template<class A>
      void transform( A const& _atoms, VariableSet &_vars,
                      Crystal::TStructure<std::string> const &_structure );

    Representation :: Representation(Crystal::TStructure<std::string> const &_structure,
                                     size_t _natoms )
    {
      // now loops over bases.
      typedef Bases< Crystal::TStructure<std::string> > t_Bases;
      t_Bases bases( _structure );
      t_Bases::const_iterator i_basis = bases.begin();
      t_Bases::const_iterator i_basis_end = bases.end();

      // first neighbor containor.
      Crystal :: Neighbors neighbors( _natoms );
      neighbors.origin = i_basis->origin + atat::rVector3d(1,0,0);
      
      // sorting container.
      std::vector<Crystal::Neighbor const*> atoms(_natoms, NULL);
      
      for(; i_basis != i_basis_end; ++i_basis )
      {
        if( not Fuzzy::is_zero( atat::norm2(i_basis->origin - neighbors.origin) ) )
        {
          neighbors.origin = i_basis->origin;
          std::vector<Crystal::Neighbor const*>::iterator i_atom = atoms.begin();
          std::vector<Crystal::Neighbor const*>::iterator i_atom_end = atoms.end();
          Crystal::Neighbors::const_iterator i_neigh = neighbors.begin( _structure );
          for(; i_atom != i_atom_end; ++i_atom, ++i_neigh ) *i_atom = &(*i_neigh);
        };
        // sorst according to basis.
        std::sort
        ( 
          atoms.begin(),
          atoms.end(),
          basis_sort(*i_basis) 
        );

        VariableSet variable_set;
        variable_set.weight = i_basis->weight;
        variable_set.variables.reserve(_natoms);
        transform( atoms, variable_set, *i_basis, _structure );
        add_(variable_set);
      };
    }

    void Representation::add_( VariableSet const &_rep )
    {
      t_Sets::iterator i_found = std::find( sets_.begin(), sets_.end(), _rep );
      if( i_found == sets_.end() ) sets_.push_back(_rep);
      else i_found->weight += _rep.weight;
    }

    template<class A>
      void transform( A const& _atoms, VariableSet &_vars, Basis const &_basis,
                      Crystal::TStructure<std::string> const &_structure )
      {
        foreach( A::pointer const ptr_atom, _atoms )
        {
          Crystal::Structure::t_Atom const atom = structure.atoms[ ptr_atom->index ];
          specie_type type(0);
          foreach( std::string const& str, structure.lattice->sites[atom.site].type )
          {
            if( str == atom.type ) break;
            ++i;
          }
          atat::rVector3d vec(atom.pos-_basis.origin); 
          switch( _vars.size() )
          {
            default: _vars.push_back( VariableSet::t_Variable(vec * _basis.z, type) );
            case 2:  _vars.push_back( VariableSet::t_Variable(vec * _basis.y, type) );
            case 1:  _vars.push_back( VariableSet::t_Variable(vec * _basis.z, type) ); break;
            case 0:  _vars.push_back( VariableSet::t_Variable(0, type) ); break;
          }
        }
      }
     
  } // namespace atomic_potential
} // namespace LaDa

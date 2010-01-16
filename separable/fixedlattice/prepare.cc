//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>


#include "prepare.h"

#include <crystal/confsplit.h>
#include <crystal/neighbors.h>
#include <atomic_potentials/bases.h>
#include <crystal/fill_structure.h>
#include <math/fuzzy.h>

namespace LaDa
{
  namespace CE
  {
    // Forward declarations.
    namespace details
    {
      void cubic_basis( types::t_unsigned _n, const Eigen::Matrix3d &_cell,
                        std::vector< Eigen::Vector3d >& _positions );
      void supercell_basis( types::t_unsigned _n, const Eigen::Matrix3d &_cell,
                            std::vector< Eigen::Vector3d >& _positions );
      void convcell_basis( types::t_unsigned _n,
                           std::vector< Eigen::Vector3d >& _positions );
    }
    
    void PosToConfs :: init_syms( Crystal::Lattice &_lat )
    {
      types::t_int N( _lat.space_group.size() );
      __ASSERT( N <= 0, "Lattice does not have symmetry operations.\n" )
      syms.reserve( N );
      for( types::t_int i(0); i < N; ++i )
      {
        if( not math::eq( atat::norm2( _lat.space_group[i].trans ), 0e0 ) ) continue;
        Eigen::Matrix3d &op = _lat.space_group[i].op;
  //     if( not math::eq( atat::det( op ), 1e0 ) ) continue;
        syms.push_back( op );
      }
    }

    void PosToConfs :: create_positions( const std::string &_bdesc )
    {
      if( _bdesc.empty() ) return;
      const boost::regex re1("(\\d+)" );
      const boost::regex re2("(\\d+)(?:\\s+)?x(?:\\s+)?"
                             "(\\d+)(?:\\s+)?x(?:\\s+)?(\\d+)" );
      const boost::regex re3("\\((\\d+)\\,(\\d+)\\)" );
      boost::match_results<std::string::const_iterator> what;

      if( boost::regex_search( _bdesc, what, re2 ) )
      {
        n.first = n.second = 0;
        Eigen::Matrix3d cell;
        cell.set_diagonal( boost::lexical_cast<types::t_real>(what.str(1)),
                           boost::lexical_cast<types::t_real>(what.str(2)),
                           boost::lexical_cast<types::t_real>(what.str(3)) );
        details::supercell_basis( 1, cell, positions );
        return;
      }
      if( boost::regex_search( _bdesc, what, re3 ) )
      {
        __TRYCODE( n.first = boost::lexical_cast< size_t >( what.str(1) );,
                   "Could not parse " << _bdesc << " -> " << what.str(1) << "\n" )
        __TRYCODE( n.second = boost::lexical_cast< size_t >( what.str(2) );,
                   "Could not parse " << _bdesc << " -> " << what.str(2) << "\n" )
        return;
      }
      if( boost::regex_search( _bdesc, what, re1 ) )
      {
        n.first = 0;
        __TRYCODE( n.second = boost::lexical_cast< size_t >( what.str(1) );,
                   "Could not parse " << _bdesc << " -> " << what.str(1) << "\n" )
        positions.clear();
        return;
      }
      __THROW_ERROR( "Could not parse --basis input: " << _bdesc << "\n" )
    }

    // Strict weak ordering functor from knowledge of basis.
    struct basis_sort
    {
      atomic_potential::Basis const& basis;
      basis_sort( atomic_potential::Basis const &_b ) : basis(_b) {}
      basis_sort( basis_sort const &_b ) : basis(_b.basis) {}
      bool operator()( Crystal::Neighbor const * const _a,
                       Crystal::Neighbor const * const _b ) const
      {
        if( math::neq( _a->distance, _b->distance ) ) return _a->distance < _b->distance;
        const types::t_real ax = _a->pos * basis.x;
        const types::t_real bx = _b->pos * basis.x;
        if( math::neq( ax, bx ) ) return ax > bx;
        const types::t_real ay = _a->pos * basis.y;
        const types::t_real by = _b->pos * basis.y;
        if( math::neq( ay, by ) ) return ay > by;
        const types::t_real az = _a->pos * basis.z;
        const types::t_real bz = _b->pos * basis.z;
        if( math::neq( az, bz ) ) return az > bz;
        return false;
      }
    }; 

    // Converts neighbor to type index.
    struct type_to_bit
    {
      Crystal::Structure const &structure;
      type_to_bit( Crystal::Structure const& _str ) : structure(_str) {}
      type_to_bit( type_to_bit const& _c ) : structure(_c.structure) {}
      size_t operator()( Crystal::Neighbor const * const _a ) const
      {
        Crystal::Structure::t_Atom const atom = structure.atoms[ _a->index ];
        return structure.lattice->convert_real_to_type_index( atom.site, atom.type );
      }
    };

    void PosToConfs :: nsites( const Crystal :: Structure &_structure,
                               t_Configurations& _confs ) const
    {
      const size_t N( _confs.size() );

      // now loops over bases.
      typedef atomic_potential::Bases< Crystal::Structure > t_Bases;
      t_Bases bases( _structure );
      t_Bases::const_iterator i_basis = bases.begin();
      t_Bases::const_iterator i_basis_end = bases.end();

      // first neighbor containor.
      Crystal :: Neighbors neighbors( n.second );
      neighbors.origin = i_basis->origin + Eigen::Vector3d(1,0,0);
      size_t index(1);
      
      // sorting container.
      std::vector<Crystal::Neighbor const*> atoms;
      size_t origin_type(0);
      
      for(; i_basis != i_basis_end; ++i_basis )
      {
        if( index != i_basis->index ) 
        {
          index = i_basis->index;

          // finds new nearest neighbors.
          neighbors.origin = i_basis->origin;
          Crystal::Neighbors::const_iterator i_neigh = neighbors.begin( _structure );
          Crystal::Neighbors::const_iterator const i_neigh_end = neighbors.end();
          atoms.clear();
          atoms.reserve(neighbors.size());
          LADA_ASSERT(neighbors.size() >= n.second, "not enough neighbors.\n")
          for(; i_neigh != i_neigh_end; ++i_neigh) atoms.push_back(&(*i_neigh));

          // finds type index at origin.
          types::t_int const sindex( _structure.atoms[index].site );
          Crystal::Structure::t_Atom::t_Type const &type = _structure.atoms[index].type;
          origin_type = _structure.lattice->convert_real_to_type_index( sindex, type );
        };
        // sorst according to basis.
        std::partial_sort
        ( 
          atoms.begin(),
          atoms.begin() + n.second - n.first,
          atoms.end(),
          basis_sort(*i_basis) 
        );
        t_Configurations :: value_type bitset;
        bitset.second = i_basis->weight;
        bitset.first.reserve(n.second - n.first + 1);
        // then copies atomic types
        bitset.first.push_back(origin_type);
        LADA_ASSERT(atoms.size() >= n.second - n.first, "not enough neighbors.\n")
        std::transform
        ( 
          atoms.begin(),
          atoms.begin() + n.second -n.first,
          std::back_inserter(bitset.first),
          type_to_bit(_structure) 
        );
        // finally adds to configurations.
        t_Configurations :: iterator i_found = _confs.begin() + N;
        t_Configurations :: iterator i_conf_end = _confs.end();
        for(; i_found != i_conf_end; ++i_found )
          if( bitset.first == i_found->first ) break;
        if( i_found == i_conf_end ) _confs.push_back( bitset );
        else i_found->second += bitset.second;
      };
    }
     
    void PosToConfs :: posbasis( const Crystal :: Structure &_structure,
                                 t_Configurations& _confs ) const
    {
      __TRYBEGIN

      using namespace boost::lambda;
      // Allocates memory.
      _confs.reserve( _confs.size() + syms.size() );
      
      types::t_real weight( 1e0 / ( syms.size() * _structure.atoms.size() ) );
      t_Positions shifted_fracs( _structure.atoms.size() );
      t_Positions fractionals( _structure.atoms.size() );
      std::transform
      ( 
        _structure.atoms.begin(), _structure.atoms.end(), fractionals.begin(),
             ret<Eigen::Matrix3d>(constant( !(_structure.cell) ))
           * bind(&Crystal::Structure::t_Atom::pos, _1) 
      );
    
      // Loops over shifts.
      typedef std::vector<Eigen::Vector3d> :: const_iterator t_shift_iterator;
      t_shift_iterator i_shift( fractionals.begin() );
      t_shift_iterator i_shift_end( fractionals.end() );
      for(; i_shift != i_shift_end; ++i_shift )
      {
        // Loops over symmetries.
        t_SymOps :: const_iterator i_op( syms.begin() );
        t_SymOps :: const_iterator i_op_end( syms.end() );
        for(; i_op != i_op_end; ++i_op )
        {
          // constructs a work array of shifted atomic positions in
          // fractional coordinates (rotation does not alter fractional coordinates).
          std::transform
          ( 
            fractionals.begin(), fractionals.end(), shifted_fracs.begin(), 
               _1 - constant( *i_shift )
          );
    
          Eigen::Matrix3d inv =  (!_structure.cell) * (*i_op);
    
          // For each basis position, finds closest atomic-position modulo
          // structure-periodicity.
          t_Bitset bitset( positions.size() );
          t_Positions :: const_iterator i_pos( positions.begin() );
          t_Positions :: const_iterator i_pos_end( positions.end() );
          for(types::t_int i=0; i_pos != i_pos_end; ++i_pos, ++i )
          {
            Eigen::Vector3d pos = inv * (*i_pos);
            t_Positions::const_iterator i_found( shifted_fracs.begin() );
            t_Positions::const_iterator i_end( shifted_fracs.end() );
            types::t_int j(0);
            for(; i_found != i_end; ++i_found, ++j )
            {
              Eigen::Vector3d a = pos - (*i_found);
              a[0] += 0.05e0; a[0] -= std::floor(a[0]); a[0] -= 0.05e0;
              a[1] += 0.05e0; a[1] -= std::floor(a[1]); a[1] -= 0.05e0;
              a[2] += 0.05e0; a[2] -= std::floor(a[2]); a[2] -= 0.05e0;
              if( math::eq( atat::norm2( a ), 0e0 ) ) break; 
            }
            __DOASSERT( i_found == shifted_fracs.end(),
                        "Could not find equivalent position.\n" ) 
            // Found the position in atomic structure.
            bitset[ i ] = Crystal::Structure::lattice
                             ->convert_real_to_type_index( 0, _structure.atoms[j].type );
          }
          // adds to configurations if similar bitset cannot be found.
          // otherwise increase weight of similar structure.
          if( _confs.size() == 0 ) _confs.push_back( t_CoefBitset( bitset, weight ) );
          else
          {
            t_Configurations :: iterator i_found;
            i_found = std::find_if
                      (
                        _confs.begin(), _confs.end(), 
                        bind( &t_CoefBitset::first, _1 ) == constant( bitset )
                      );
            if( i_found == _confs.end() ) 
              _confs.push_back( t_CoefBitset( bitset, weight ) );
            else i_found->second += weight;
          }
        } // loop over symmetry operations.
    
      } // loop over origin-shifts 
      __TRYEND(,"Error while creating configurations.\n" )
    }
    
    namespace details
    {
      void cubic_basis( types::t_unsigned _n, const Eigen::Matrix3d &_cell,
                        std::vector< Eigen::Vector3d >& _positions )
      {
        _positions.clear();
        _positions.reserve( _n*_n*_n );
        for( types::t_unsigned i = 0; i < _n; ++i )
          for( types::t_unsigned j = 0; j < _n; ++j )
            for( types::t_unsigned k = 0; k < _n; ++k )
            {
              Eigen::Vector3d pos( i, j, k );
              _positions.push_back( _cell * pos );
            }
        namespace bl = boost::lambda;
        std::sort
        ( 
          _positions.begin(), _positions.end(),
           bl::_1 * bl::_1  > bl::_2 * bl::_2
        );
      }
      void supercell_basis( types::t_unsigned _n, const Eigen::Matrix3d &_cell,
                            std::vector< Eigen::Vector3d >& _positions )
      {
        __DEBUGTRYBEGIN
        namespace bl = boost::lambda;

        Crystal::Structure structure;
        Eigen::Matrix3d mat;
        mat.identity(); mat =  mat * (types::t_real)_n;
        structure.cell = _cell * mat;
        __ASSERT( not Crystal::Structure::lattice,
                  "Lattice of structure has not beens set.\n" )
        Crystal :: fill_structure( structure );
        std::for_each
        (
          structure.atoms.begin(), structure.atoms.end(),
          ( 
            bl::bind( &Crystal::Structure::t_Atom::type, bl::_1 ) = 1e0,
            bl::bind( &Crystal::Structure::t_Atom::site, bl::_1 ) = 0
          )
        );

        _positions.clear();
        _positions.resize( structure.atoms.size() );
        std::transform
        (
          structure.atoms.begin(), structure.atoms.end(), _positions.begin(),
          bl::bind( &Crystal::Structure::t_Atom::pos, bl::_1 )
        );
        std::sort
        ( 
          _positions.begin(), _positions.end(),
           bl::_1 * bl::_1  > bl::_2 * bl::_2
        );
        __DEBUGTRYEND(, "Error while creating super-cell basis.\n" )
      }
      void convcell_basis( types::t_unsigned _n, 
                           std::vector< Eigen::Vector3d >& _positions )
      {
        __DEBUGTRYBEGIN
        __ASSERT( not Crystal::Structure::lattice, 
                  "Lattice has not been set.\n" )
        Eigen::Matrix3d mult;
        // assume fcc
        if( math::eq( Crystal::Structure::lattice->cell.x[0][0], 0e0 ) ) 
        {
          mult.x[0][0] = -1e0; mult.x[1][0] =  1e0; mult.x[2][0] =  1e0; 
          mult.x[0][1] =  1e0; mult.x[1][1] = -1e0; mult.x[2][1] =  1e0; 
          mult.x[0][2] =  1e0; mult.x[1][2] =  1e0; mult.x[2][2] = -1e0; 
        }
        else // assume bcc
        {
          mult.x[0][0] = 0e0; mult.x[1][0] = 1e0; mult.x[2][0] = 1e0; 
          mult.x[0][1] = 1e0; mult.x[1][1] = 0e0; mult.x[2][1] = 1e0; 
          mult.x[0][2] = 1e0; mult.x[1][2] = 1e0; mult.x[2][2] = 0e0; 
        }
        mult = Crystal::Structure::lattice->cell * mult;
        supercell_basis( _n, mult, _positions );
        __DEBUGTRYEND(, "Error while creating conventional-cell basis.\n" )
      }
    }

  }
} // namespace LaDa

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
#include <crystal/fill_structure.h>

namespace LaDa
{
  namespace CE
  {
    // Forward declarations.
    namespace details
    {
      void cubic_basis( types::t_unsigned _n, const atat::rMatrix3d &_cell,
                        std::vector< atat::rVector3d >& _positions );
      void supercell_basis( types::t_unsigned _n, const atat::rMatrix3d &_cell,
                            std::vector< atat::rVector3d >& _positions );
      void convcell_basis( types::t_unsigned _n,
                           std::vector< atat::rVector3d >& _positions );
    }
    
    void PosToConfs :: init_syms( Crystal::Lattice &_lat )
    {
      types::t_int N( _lat.space_group.point_op.getSize() );
      __ASSERT( N <= 0, "Lattice does not have symmetry operations.\n" )
      syms.reserve( N );
      for( types::t_int i(0); i < N; ++i )
      {
        if( not Fuzzy::eq( atat::norm2( _lat.space_group.trans(i) ), 0e0 ) ) continue;
        atat::rMatrix3d &op = _lat.space_group.point_op(i);
  //     if( not Fuzzy::eq( atat::det( op ), 1e0 ) ) continue;
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
        atat::rMatrix3d cell;
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

//   void PosToConfs :: nsites( const Crystal :: Structure &_structure,
//                              t_Configurations& _confs ) const
//   {
//     Crystal::SplitIntoConfs split;
//     split( _structure, n );
//     const size_t N( _confs.size() );
//     _confs.reserve( N + split.configurations().size() );
//     typedef Crystal::SplitIntoConfs::t_Configurations tt_confs;
//     typedef tt_confs :: value_type tt_conf;
//     typedef tt_conf :: first_type tt_bits;
//     foreach( const tt_conf & _atomconf,  split.configurations() )
//     {
//       t_Configurations :: value_type bitset;
//       bitset.second = _atomconf.second;
//       bitset.first.resize( _atomconf.first.size() );
//       typedef t_Configurations :: value_type :: first_type t_bits;
//       t_bits :: iterator i_bit = bitset.first.begin();
//       foreach( const tt_bits :: value_type &bit, _atomconf.first )
//       {
//         const Crystal::Structure::t_Atom &atom( _structure.atoms[ bit.second ] );
//         const Crystal::Structure::t_Atom::t_Type type( atom.type );
//         *i_bit = Crystal::Structure::lattice
//                      ->convert_real_to_type_index( atom.site, type );
//         ++i_bit;
//       }
//       // adds to configurations.
//       t_Configurations :: iterator i_found = _confs.begin() + N;
//       t_Configurations :: iterator i_conf_end = _confs.end();
//       for(; i_found != i_conf_end; ++i_found )
//         if( bitset.first == i_found->first ) break;
//       if( i_found == i_conf_end ) _confs.push_back( bitset );
//       else i_found->second += bitset.second;
//     }
//     split.clear();
//   }

    void PosToConfs :: nsites( const Crystal :: Structure &_structure,
                               t_Configurations& _confs ) const
    {
      struct basis_sort
      {
        configuration::basis const& basis;
        basis_sort( configuration::basis const &_b ) : basis(_b) {}
        basis_sort( basis_sort const &_b ) : basis(_b.basis) {}
        bool operator()( atat::rVector3d const *_a, atat::rVector3d const *_b ) 
      }; 
      // Copy atomic positions as pointers.
      std::vector< Crystal::Structure::t_Atom*> atoms;
      atoms.reserve(_structure.atoms.size());
      foreach( Crystal::Structure::t_Atom const & atom, _structure.atoms )
        atoms.push_back( &atom );

      // now loops over bases.
      CE::configurations::Bases bases( _structure );
      CE::configurations::Bases::const_iterator i_basis = bases.begin();
      CE::configurations::Bases::const_iterator i_basis_end = bases.end();
      for(; i_basis != i_basis_end; ++i_basis )
      {
        // sorst according to basis.
        std::partial_sort( atoms.begin(), atoms.begin() + n.second, basis_sort(*i_basis) );
        t_Configurations :: value_type bitset;
        bitset.second = i_basis->weight;
        // then copies atomic types
        std::transform
        ( 
          atoms.begin() + n.first, atoms.begin() + n.second,
          std::back_inserter(bitset),
          type_to_bit(_structure.lattice) 
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
             ret<atat::rMatrix3d>(constant( !(_structure.cell) ))
           * bind(&Crystal::Structure::t_Atom::pos, _1) 
      );
    
      // Loops over shifts.
      typedef std::vector<atat::rVector3d> :: const_iterator t_shift_iterator;
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
    
          atat::rMatrix3d inv =  (!_structure.cell) * (*i_op);
    
          // For each basis position, finds closest atomic-position modulo
          // structure-periodicity.
          t_Bitset bitset( positions.size() );
          t_Positions :: const_iterator i_pos( positions.begin() );
          t_Positions :: const_iterator i_pos_end( positions.end() );
          for(types::t_int i=0; i_pos != i_pos_end; ++i_pos, ++i )
          {
            atat::rVector3d pos = inv * (*i_pos);
            t_Positions::const_iterator i_found( shifted_fracs.begin() );
            t_Positions::const_iterator i_end( shifted_fracs.end() );
            types::t_int j(0);
            for(; i_found != i_end; ++i_found, ++j )
            {
              atat::rVector3d a = pos - (*i_found);
              a[0] += 0.05e0; a[0] -= std::floor(a[0]); a[0] -= 0.05e0;
              a[1] += 0.05e0; a[1] -= std::floor(a[1]); a[1] -= 0.05e0;
              a[2] += 0.05e0; a[2] -= std::floor(a[2]); a[2] -= 0.05e0;
              if( Fuzzy::eq( atat::norm2( a ), 0e0 ) ) break; 
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
      void cubic_basis( types::t_unsigned _n, const atat::rMatrix3d &_cell,
                        std::vector< atat::rVector3d >& _positions )
      {
        _positions.clear();
        _positions.reserve( _n*_n*_n );
        for( types::t_unsigned i = 0; i < _n; ++i )
          for( types::t_unsigned j = 0; j < _n; ++j )
            for( types::t_unsigned k = 0; k < _n; ++k )
            {
              atat::rVector3d pos( i, j, k );
              _positions.push_back( _cell * pos );
            }
        namespace bl = boost::lambda;
        std::sort
        ( 
          _positions.begin(), _positions.end(),
           bl::_1 * bl::_1  > bl::_2 * bl::_2
        );
      }
      void supercell_basis( types::t_unsigned _n, const atat::rMatrix3d &_cell,
                            std::vector< atat::rVector3d >& _positions )
      {
        __DEBUGTRYBEGIN
        namespace bl = boost::lambda;

        Crystal::Structure structure;
        atat::rMatrix3d mat;
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
                           std::vector< atat::rVector3d >& _positions )
      {
        __DEBUGTRYBEGIN
        __ASSERT( not Crystal::Structure::lattice, 
                  "Lattice has not been set.\n" )
        atat::rMatrix3d mult;
        // assume fcc
        if( Fuzzy::eq( Crystal::Structure::lattice->cell.x[0][0], 0e0 ) ) 
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

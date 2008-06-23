//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "ce.h"

namespace CE
{
  Separables :: Separables   ( types::t_unsigned _rank,
                               types::t_unsigned _size,
                               std::string _type )
                           : basis_size( _size ), basis_type( _type )
  {
    set_rank( _rank );
    set_basis( _size, _type );
  }

  void Separables :: set_basis( types::t_unsigned _size, std::string _type )
  {
    __ASSERT( _type.compare("cube"), "Unknown basis type.\n" )
    __ASSERT( Crystal::Structure::lattice == NULL, "Lattice type has not been set.\n" )
    basis_size = _size;
    basis_type = _type;
    details::cubic_basis( basis_size, Crystal::Structure::lattice->cell, positions );
    t_Basis :: iterator i_sep = basis.begin();
    t_Basis :: iterator i_sep_end = basis.begin();
    for(; i_sep != i_sep_end; ++i_sep )
    {
      i_sep->basis.resize( positions.size() );
      i_sep->coefs.resize( positions.size(), 0 );
    }
  }



  void SymSeparables :: init_syms( Crystal::Lattice &_lat )
  {
    types::t_int N( _lat.space_group.point_op.getSize() );
    syms.reserve( N );
    for( types::t_int i(0); i < N; ++i )
    {
      if( not Fuzzy::eq( atat::norm2( _lat.space_group.trans(i) ), 0e0 ) ) continue;
      atat::rMatrix3d &op = _lat.space_group.point_op(i);
      if( not Fuzzy::eq( atat::det( op ), 1e0 ) ) continue;
      syms.push_back( op );
    }
  }

  SymSeparables :: t_Configurations* 
    SymSeparables :: configurations( Crystal :: Structure &_structure )
    {
      t_Configurations *result(NULL);
      try
      {
        // Allocates memory.
        result = new t_Configurations;
        result->reserve( syms.size() );
        
        types::t_real weight( 1e0 / ( syms.size() * _structure.atoms.size() ) );
        work.resize( _structure.atoms.size() );
        atat::rMatrix3d inv = !_structure.cell;
        // Loops over shifts.
        typedef Crystal::Structure::t_Atoms::const_iterator t_shift_iterator;
        t_shift_iterator i_shift( _structure.atoms.begin() );
        t_shift_iterator i_shift_end( _structure.atoms.end() );
        for(; i_shift != i_shift_end; ++i_shift )
        {
          // Loops over symmetries.
          t_SymOps :: const_iterator i_op( syms.begin() );
          t_SymOps :: const_iterator i_op_end( syms.end() );
          for(; i_op != i_op_end; ++i_op )
          {
            // constructs a work array of shifted, rotated atomic positions in
            // fractional coordinates.
            atat::rVector3d shift = inv * i_shift->pos;
            atat::rMatrix3d op = inv * (*i_op );
            using namespace boost::lambda;
            std::transform
            ( 
              _structure.atoms.begin(), _structure.atoms.end(), work.begin(),
              ret<atat::rVector3d>(
                ret<atat::rVector3d>( constant(op) * bind(&Crystal::Structure::t_Atom::pos, _1) )
                 -  constant( shift ) )
            );

            // For each basis position, finds closest atomic-position modulo
            // structure-periodicity.
            t_Bitset bitset( basis.size() );
            t_Basis :: const_iterator i_pos( basis.begin() );
            t_Basis :: const_iterator i_pos_end( basis.end() );
            for(types::t_int i=0; i_pos != i_pos_end; ++i_pos, ++i )
            {
              atat::rVector3d pos = inv * (*i_pos);
              t_Basis::const_iterator i_found( work.begin() );
              t_Basis::const_iterator i_end( work.end() );
              types::t_int j(0);
              for(; i_found != i_end; ++i_found, ++j )
              {
                atat::rVector3d a = pos - (*i_found);
                a[0] += 0.5e0; a[0] -= std::floor(a[0]); a[0] -= 0.5e0;
                a[1] += 0.5e0; a[1] -= std::floor(a[1]); a[1] -= 0.5e0;
                a[2] += 0.5e0; a[2] -= std::floor(a[2]); a[2] -= 0.5e0;
                if( Fuzzy::eq( atat::norm2( a ), 0e0 ) ) break; 
              }
              __DOASSERT( i_found == work.end(), "Could not find equivalent position.\n" ) 
              // Found the position in atomic structure.
              bitset[ i ] = Fuzzy::eq( _structure.atoms[j].type, -1e0 );
            }
            // adds to configurations if similar bitset cannot be found.
            // otherwise increase weight of similar structure.
            if( result->size() == 0 ) result->push_back( t_CoefBitset( bitset, weight ) );
            else
            {
              t_Configurations :: iterator i_found;
              i_found = std::find_if
                        (
                          result->begin(), result->end(), 
                          bind( &t_CoefBitset::first, _1 ) == constant( bitset )
                        );
              if( i_found == result->end() ) 
                result->push_back( t_CoefBitset( bitset, weight ) );
              else i_found->second += weight;
            }
          } // loop over symmetry operations.

        } // loop over origin-shifts 

      }  // try
      __CATCHCODE( if( result ) delete result;,
                   "Error while creating configurations.\n" )
    
    }
  
    types::t_real SymSeparables :: operator()( t_Configurations &_configs, 
                                               const Separables &_func ) const
    {
      using namespace boost::lambda;
      types::t_real result(0);
      std::for_each
      (
        _configs.begin(), _configs.end(),
        var(result) +=
          ret<t_CoefBitset::second_type>
             (
                 bind( &t_CoefBitset::second, _1 )
               * bind<t_CoefBitset::second_type>( _func, bind( &t_CoefBitset::first, _1 ) )
             )
      );
      return result; 
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
    }
  }

}

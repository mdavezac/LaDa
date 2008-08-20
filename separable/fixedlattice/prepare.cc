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
      if( not Fuzzy::eq( atat::det( op ), 1e0 ) ) continue;
      syms.push_back( op );
    }
  }

  void PosToConfs :: create_positions( const std::string &_bdesc )
  {
    if( _bdesc.empty() ) return;
    const boost::regex re("(\\d+)(?:\\s+)?x(?:\\s+)?"
                          "(\\d+)(?:\\s+)?x(?:\\s+)?(\\d+)" );
    boost::match_results<std::string::const_iterator> what;
    __DOASSERT( not boost::regex_search( _bdesc, what, re ),
                "Could not parse --basis input: " << _bdesc << "\n" )
    atat::rMatrix3d cell;
    cell.set_diagonal( boost::lexical_cast<types::t_real>(what.str(1)),
                       boost::lexical_cast<types::t_real>(what.str(2)),
                       boost::lexical_cast<types::t_real>(what.str(3)) );
    details::supercell_basis( 1, cell, positions );
  }

   
  void PosToConfs :: operator()( const Crystal :: Structure &_structure,
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
           ret<atat::rMatrix3d>(constant( !(~_structure.cell) ))
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
  
        atat::rMatrix3d inv = !( (*i_op) * (~_structure.cell)  );
  
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
          bitset[ i ] = Fuzzy::eq( _structure.atoms[j].type, -1e0 );
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

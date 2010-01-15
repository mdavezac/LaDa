//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <Eigen/LU>

#include <crystal/fill_structure.h>

#include "ce.h"

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
    
    Separables :: Separables   ( types::t_unsigned _rank,
                                 types::t_unsigned _size,
                                 const std::string &_type )
                             : t_Base(), basis_size( _size ), basis_type( _type ), 
                               name("Sum of Separables")
    {
      set_rank( _rank );
      set_basis( _size, _type );
    }
    Separables :: Separables   ( types::t_unsigned _rank,
                                 const Eigen::Matrix3d &_cell,
                                 const std::string &_type )
                             : t_Base(), basis_type( _type ), 
                               name("Sum of Separables")
    {
      set_rank( _rank );
      set_basis( _cell, _type );
    }

    void Separables :: set_basis( const Eigen::Matrix3d &_cell, const std::string &_type )
    {
      __ASSERT( Crystal::Structure::lattice == NULL,
                "Lattice type has not been set.\n" )
      details::supercell_basis( 1, _cell, positions );
      basis_type = _type;
      basis_size = positions.size();
      rename();
    }
    void Separables :: set_basis( types::t_unsigned _size, const std::string &_type )
    {
      __ASSERT( _type.compare("cube") and _type.compare("conv"),
                "Unknown basis type " << _type << ".\n")
      __ASSERT( Crystal::Structure::lattice == NULL,
                "Lattice type has not been set.\n" )
      if( basis_type.compare("conv") == 0 ) 
        details::convcell_basis( _size, positions );
      else details::cubic_basis( _size,
                                 Crystal::Structure::lattice->cell,
                                 positions );
      basis_type = _type;
      basis_size = positions.size();
      rename();
    }
    void Separables :: set_rank( types::t_unsigned _rank ) 
    {
      basis.resize( _rank );
      coefs.resize( _rank, 1e0 );
      namespace bl = boost::lambda;
      std::for_each // loop over ranks.
      ( 
        basis.begin(), basis.end(),
        bl::bind
        (
          &t_Basis::value_type::set_name,
          bl::_1, bl::constant("Rank Function") 
        ) // end of bl::bind
      ); // end of loop over ranks.
    }

    void Separables::rename() 
    {
      // loop over ranks.
      t_Basis :: iterator i_sep = basis.begin();
      t_Basis :: iterator i_sep_end = basis.end();
      for(; i_sep != i_sep_end; ++i_sep )
      {
        i_sep->basis.resize( positions.size() );
        i_sep->coefs.resize( positions.size(), 0 );
        namespace bl = boost::lambda;
        std::for_each // loop over dimensions
        ( 
          i_sep->basis.begin(), i_sep->basis.end(),
          bl::bind
          (
            &t_Basis::value_type::t_Basis::value_type::set_name,
            bl::_1, bl::constant("Dimensional Function")
          ) // end of bl::bind
        ); // end of loop over dimensions
      } // end of loop over ranks
    }



    SymSeparables :: SymSeparables ( Separables &_sep ) : basis( _sep.positions )
    {
      if ( not Crystal::Structure::lattice ) return;
      
      if( not Crystal::Structure::lattice->space_group.size() )
        Crystal::Structure::lattice->find_space_group();
      try{ init_syms( *Crystal::Structure::lattice ); }
      catch(...){}
    }

    void SymSeparables :: init_syms( Crystal::Lattice &_lat )
    {
      types::t_int N( _lat.space_group.size() );
      LADA_ASSERT( N >= 0, "Lattice does not have symmetry operations.\n" )
      syms.reserve( N );
      for( types::t_int i(0); i < N; ++i )
      {
        if( not math::is_zero(_lat.space_group[i].trans.squaredNorm()) ) continue;
        Eigen::Matrix3d &op = _lat.space_group[i].op;
        if( not math::is_zero(op.determinant() - 1e0) ) continue;
        syms.push_back( op );
      }
    }

     
    void SymSeparables :: configurations( const Crystal :: Structure &_structure,
                                          SymSeparables :: t_Configurations& _confs ) const
    {
      try
      {
        using namespace boost::lambda;
        // Allocates memory.
        _confs.reserve( _confs.size() + syms.size() );
        
        types::t_real weight( 1e0 / ( syms.size() * _structure.atoms.size() ) );
        t_Basis shifted_fracs( _structure.atoms.size() );
        t_Basis fractionals( _structure.atoms.size() );
        std::transform
        ( 
          _structure.atoms.begin(), _structure.atoms.end(), fractionals.begin(),
               ret<Eigen::Matrix3d>(constant( !(~_structure.cell) ))
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

            Eigen::Matrix3d inv = !( (*i_op) * (~_structure.cell)  );

            // For each basis position, finds closest atomic-position modulo
            // structure-periodicity.
            t_Bitset bitset( basis.size() );
            t_Basis :: const_iterator i_pos( basis.begin() );
            t_Basis :: const_iterator i_pos_end( basis.end() );
            for(types::t_int i=0; i_pos != i_pos_end; ++i_pos, ++i )
            {
              Eigen::Vector3d pos = inv * (*i_pos);
              t_Basis::const_iterator i_found( shifted_fracs.begin() );
              t_Basis::const_iterator i_end( shifted_fracs.end() );
              types::t_int j(0);
              for(; i_found != i_end; ++i_found, ++j )
              {
                Eigen::Vector3d a = pos - (*i_found);
                a[0] += 0.05e0; a[0] -= std::floor(a[0]); a[0] -= 0.05e0;
                a[1] += 0.05e0; a[1] -= std::floor(a[1]); a[1] -= 0.05e0;
                a[2] += 0.05e0; a[2] -= std::floor(a[2]); a[2] -= 0.05e0;
                if( math::is_zero( a.squaredNorm(), 0e0 ) ) break; 
              }
              __DOASSERT( i_found == shifted_fracs.end(), "Could not find equivalent position.\n" ) 
              // Found the position in atomic structure.
              bitset[ i ] = math::eq( _structure.atoms[j].type, -1e0 );
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
      }  // try
      __CATCHCODE(,"Error while creating configurations.\n" )
    
    }
    
    types::t_real SymSeparables :: operator()( t_Configurations &_configs, 
                                               const Separables &_func ) const
    {
      namespace bl = boost::lambda;
      types::t_real result(0);
      t_Configurations :: const_iterator i_conf = _configs.begin();
      t_Configurations :: const_iterator i_conf_end = _configs.end();
      for(; i_conf != i_conf_end; ++i_conf )
        result += i_conf->second * _func( i_conf->first.begin(),
                                          i_conf->first.end()   );
      return result; 
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

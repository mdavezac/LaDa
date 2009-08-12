//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_READ_POSCAR_H_
#define _LADA_CRYSTAL_READ_POSCAR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <ostream>
#include <fstream>
#include <string>
#include <numeric>

#include <boost/filesystem/operations.hpp>
#include <boost/spirit/include/classic_numerics.hpp>
#include <boost/spirit/include/classic_primitives.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/classic_operators.hpp>
#include <boost/spirit/include/classic_kleene_star.hpp>
#include <boost/spirit/include/classic_actions.hpp>

#include <print/manip.h>
#include <opt/types.h>

#include "structure.h"

namespace LaDa 
{
  namespace Crystal
  {

    //! Reads structure in NREL format.
    template<class T_TYPE> 
      void read_poscar( TStructure<T_TYPE> &_structure, 
                        const boost::filesystem::path &_path,
                        const std::vector< T_TYPE >& _types );

    template<class T_TYPE> 
      void read_poscar( TStructure<T_TYPE> &_structure, 
                        const boost::filesystem::path &_path,
                        const std::vector< T_TYPE >& _types )
      {
        __TRYBEGIN
        
        namespace bsc = boost::spirit::classic;
        namespace fs = boost::filesystem;  
        __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
        __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                    _path << " is neither a regulare file nor a system link.\n" )
        std::ifstream file( _path.string().c_str(), std::ifstream::in );
        std::string line;
        __DOASSERT( file.bad(), "Could not open file " + _path.string() + "\n" );
        
        // name
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        _structure.name = Print::StripEdges( line );
        // scale
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
            bsc::ureal_p[ bsc::assign_a( _structure.scale ) ],
            bsc::space_p
          ).hit,
          "Could not parse cell.\n" 
        )
        // cell.
        for( size_t i(0); i < 3; ++i )
        {
          std::getline( file, line );
          __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
          __DOASSERT
          (
            not bsc::parse
            (
              line.c_str(),
                  bsc::real_p[ bsc::assign_a( _structure.cell(i, 0) ) ] 
               >> bsc::real_p[ bsc::assign_a( _structure.cell(i, 1) ) ] 
               >> bsc::real_p[ bsc::assign_a( _structure.cell(i, 2) ) ],
               bsc::space_p
            ).hit,
            "Could not parse cell.\n" 
          )
        }
        // number of atoms
        std::vector< size_t > nbatoms;
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
            *( bsc::uint_p[ bsc::push_back_a( nbatoms ) ] ),
             bsc::space_p
          ).hit,
          "Could not parse number of atoms per specie.\n" 
        )

        // Cartesian or direct.
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        line = Print :: StripEdges(line);
        bool direct;
        if( line[0] == 'c' or line[0] == 'C' ) direct = false;
        else if( line[0] == 'd' or line[0] == 'D' ) direct = true;
        else __THROW_ERROR( "Unknown cartesian/direct switch.\n" )
        
        // Atomic positions.
        __DOASSERT( nbatoms.size() != _types.size(), 
                       "Number of species in POSCAR does not match input: " 
                    << nbatoms.size() << " " << _types.size() << ".\n" )
        const size_t nb( std::accumulate( nbatoms.begin(), nbatoms.end(), 0 ) );
        typename std::vector<T_TYPE> :: const_iterator i_type( _types.begin() );
        std::vector<size_t> :: const_iterator i_nb( nbatoms.begin() );
        typename TStructure<T_TYPE> :: t_Atom atom;
        atom.freeze = TStructure<T_TYPE> :: t_Atom :: FREEZE_NONE;
        atom.site = -1;
        for( size_t i(0), j(1); i < nb; ++i, ++j )
        {
          __ASSERT( i_nb == nbatoms.end(), "Unexpected end of vector.\n" )
          __ASSERT( i_type == _types.end(), "Unexpected end of vector.\n" )
          std::getline( file, line );
          __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
          __DOASSERT
          (
            not bsc::parse
            (
              line.c_str(),
            
                  bsc::real_p[ bsc::assign_a( atom.pos[0] ) ] 
               >> bsc::real_p[ bsc::assign_a( atom.pos[1] ) ] 
               >> bsc::real_p[ bsc::assign_a( atom.pos[2] ) ],
               bsc::space_p
            ).hit,
            "Could not parse atomic position " <<  i << ".\n" 
          )
          if( direct ) atom.pos = _structure.cell * atom.pos; 
          atom.type = *i_type;
          _structure.atoms.push_back( atom );
          if( j != *i_nb ) continue;
          ++i_nb; ++i_type; j = 0;
        }
        __TRYEND(, "Could not parse POSCAR " + _path.string() + "\n" )
      }

  } // namespace Crystal

} // namespace LaDa

#endif

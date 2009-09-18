//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/ref.hpp>

#include <boost/spirit/include/classic_primitives.hpp>
#include <boost/spirit/include/classic_numerics.hpp>
#include <boost/spirit/include/classic_parser.hpp>
#include <boost/spirit/include/classic_actions.hpp>
#include <boost/spirit/include/classic_kleene_star.hpp>
#include <boost/spirit/include/classic_optional.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_sequence.hpp>
#include <boost/spirit/include/classic_operators.hpp>

#include <opt/ndim_iterator.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>

#include <atat/misc.h>

#include "structure.h"
#include "read_structure.h"



namespace LaDa
{

  namespace Crystal 
  {
    void read_structure( Structure &_struct,
                         const boost::filesystem::path &_path )
    {
      __TRYBEGIN

      namespace fs = boost::filesystem;  
      __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                  _path << " is neither a regulare file nor a system link.\n" )
      std::ifstream file( _path.string().c_str(), std::ifstream::in );
      std::string line;

      std::getline( file, line ); // name and inconsequential data.

      _struct.name = _path.string();
      types::t_int N;  // number of atoms;
      std::getline( file, line ); // name and inconsequential data.
      std::istringstream sstr( line );
      sstr >> N;
      // cell 
      for(size_t i(0); i < 3; ++i )
      {
        __ASSERT( file.eof(),
                  "Reached unexpected end of file: " << _path << ".\n" )
        std::getline( file, line );
        sstr.str( line ); sstr.seekg (0, std::ios::beg); sstr.clear();
        sstr >> _struct.cell.x[0][i]
             >> _struct.cell.x[1][i]
             >> _struct.cell.x[2][i];
      }
      _struct.freeze = Crystal::Structure::FREEZE_NONE;
      LADA_DOASSERT( atat::is_int((!_struct.cell) * _struct.cell),
                     "Structure cell is not supercell of lattice." )
      // now atoms.
      types::t_int nfound(0);
      while( nfound < N and file.good() )
      {
        __ASSERT( file.eof(),
                  "Reached unexpected end of file: " << _path << ".\n" )
        std::getline( file, line );
        sstr.str( line ); sstr.seekg (0, std::ios::beg); sstr.clear();
        Crystal::Structure::t_Atom a;
        types::t_int type;
        sstr >> type;
        ++nfound;
        if( type != 1 and type != 2 ) continue;
        a.type = ( type == 1 ) ? -1.e0: 1.e0;
        sstr >> a.pos.x[0] >> a.pos.x[1]  >> a.pos.x[2];
        a.freeze = Structure::t_Atom::FREEZE_NONE;
        a.site = 0;
        _struct.atoms.push_back(a);
      }
      __ASSERT( nfound != N,    "Could find only " << nfound << " of " 
                             << N << " atoms in " << _path << ".\n" )
        
      __TRYEND(, "Error while reading " << _path << "\n" )
    }

    void read_ce_structures( const boost::filesystem::path &_dir,
                             std::vector<Crystal::Structure> &_structures )
    {
      __TRYBEGIN
      namespace fs = boost::filesystem;
      namespace bsc = boost::spirit::classic;

      // First finds directory of LDAs.dat.
      __DOASSERT( not fs::exists( _dir ), _dir << " does not exist.\n" );
      __DOASSERT( not ( fs::is_regular( _dir ) or fs::is_symlink( _dir ) ),
                  _dir << " is a not a valid file.\n" );
      fs::path dir( _dir.branch_path()  );

      // then starts reading file.
      std::ifstream ldas( _dir.string().c_str(), std::ifstream::in );
      std::string line;
      while( std::getline( ldas, line ) )
      {
        Crystal :: Structure structure;
        bsc::parse_info<> info = bsc::parse
        (
          line.c_str(),

          !(*bsc::space_p)
            >> ( *(
                    bsc::alnum_p | bsc::punct_p
               ) )[ bsc::assign_a( structure.name ) ]
            >> *bsc::space_p
            >> bsc::real_p[ bsc::assign_a( structure.energy ) ]
        );
        if( not info.hit ) continue;

        Crystal :: read_structure( structure, dir / structure.name );
        _structures.push_back(structure);
      }

      __TRYEND(,"Error while parsing " << _dir << " and structures.\n" )
    }

    bool read_pifile_structure( std::istream &_sstr,
                                Crystal::Structure &_structure )
    {
      __DEBUGTRYBEGIN
      // finds first line for structure.
      std::string line;
      do
      { 
        if( _sstr.eof() ) return false;
        std::getline( _sstr, line ); 
      } 
      while( line.find( "NO." ) == std::string::npos );

      // Tokenize first line.
      boost::char_separator<char> sep(" ");
      typedef boost::tokenizer< boost::char_separator<char> > t_Tokenizer;
      t_Tokenizer notoken(line, sep);
      t_Tokenizer::const_iterator i_tok = notoken.begin();
      t_Tokenizer::const_iterator i_tok_end = notoken.end();
      // read name.
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      _structure.name = *i_tok;
      // read size.
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      const size_t N( boost::lexical_cast<size_t>( *i_tok ) );
      // read decoration.
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      const types::t_int decoration( boost::lexical_cast<types::t_int>( *i_tok ) );
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
      // read cell.
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j, ++i_tok )
        {
          __DOASSERT( i_tok == i_tok_end, "Unexpected end-of-line.\n" )
          _structure.cell(j,i)
            = boost::lexical_cast<types::t_real>( *i_tok ) * 0.5e0;
        }

      // read atoms position.
      // find first line of atomic basis.
      _structure.atoms.clear();
      do
      { 
        if( _sstr.eof() ) return false;
        std::getline( _sstr, line ); 
      } 
      while( line.find( "BASIS" ) == std::string::npos );

      bool is_first = true;
      while( _structure.atoms.size() < N ) 
      {
        __DOASSERT( _sstr.eof(), "Unexpected end-of-file.\n" )
        t_Tokenizer basistoken(line, sep);
        i_tok = basistoken.begin();
        i_tok_end = basistoken.end();
        if( is_first ) 
        {
          is_first = false;
          __DOASSERT( (++i_tok) == i_tok_end, "Unexpected end-of-line.\n" )
        }
        Crystal::Structure :: t_Atom atom;
        atom.site = 0;
        while( i_tok != i_tok_end )
        {
          atom.type = ( decoration >> _structure.atoms.size() ) % 2 ? -1e0: 1e0;
          for( size_t i(0); i < 3; ++i, ++i_tok )
          {
            __DOASSERT( i_tok == i_tok_end, "Unexpected end-of-line.\n" )
            atom.pos(i) = boost::lexical_cast<Crystal::Structure::t_Atom::t_Type> 
                                             ( *i_tok ) * 0.5e0;
          }
          _structure.atoms.push_back( atom );
          if( _structure.atoms.size() == N ) break;
        }
        std::getline( _sstr, line );
      }
      __DOASSERT( _structure.atoms.size() != N, "Read too many atoms...\n" )

      _structure.scale = 1e0;
      _structure.k_vecs.clear();
      _structure.find_k_vectors();
      return true;
      __DEBUGTRYEND(, "Error while reading from pifile.\n" )
    }
  } // namespace Crystal

} // namespace LaDa

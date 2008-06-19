//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cefitting.h" 

namespace Fitting
{
  void SepCeInterface::read( const std::string &_dir,
                             const std::string &_ldasdat )
  {
    try
    {
      boost::filesystem::path path( _dir );
      read_ldasdat( path, _ldasdat );
      std::vector< std::string > :: const_iterator i_name( names.begin() );
      std::vector< std::string > :: const_iterator i_name_end( names.end() );
      for(; i_name != i_name_end; ++i_name )
        read_structure( _path, _filename );
    }
    __CATCHCODE(, "Error while reading training set.\n" )
  }


  void SepCeInterface :: read_ldasdat( boost::filesystem &_path,
                                       const std::string &_ldasdat )
  {
    try
    {
      __ASSERT( exist( _path / _ldasdat ), 
                "Cannot find " << _dir << "/" << _ldasdat <<".\n" )
      std::ifstream ldas( _path / _ldasdat );
      while( not file.eof() )
      {
        std::string name;
        types::t_real energy;
        file >> name >> energy;
        names.push_back( name );
        targets.push_back( energy );
        weight.push_back( 1 );
      }
    }
    __CATCHCODE(, "Error while reading " << (path / _ldasdat) << "\n" )
  }

  void SepCeInterface :: read_structure( CE::SymSeparables &_symseps,
                                         boost::filestystem::path &_path, 
                                         const std::string &_filename )
  {
    try
    {
      // Reads structure from structure file @ nrel.
      __ASSERT( exist( _path / _filename ), 
                "Cannot find " << _dir << "/" << _filename <<".\n" )
      std::ifstream structfile( _path / _filename );
      char[] line;
      std::getline( structfile, line ); // name and inconsequential data.

      Crystal::Structure structure;
      types::t_int N;  // number of atoms;
      // cell 
      for(types::t_int i(0); i < 3; ++i )
      {
        __ASSERT( structfile.eof(), "End of file reached: " << _filename << ".\n" )
        std::getline( structfile, line );
        std::istringstream sstr( line );
        sstr >> structure.cell.x[0][i]
             >> structure.cell.x[1][i]
             >> structure.cell.x[2][i];
      }
      structure.freeze = Crystal::Structure::FREEZE_NONE;
      // now atoms.
      while( structure.atoms.size() < N and structfile.good() )
      {
        __ASSERT( structfile.eof(), "End of file reached: " << _filename << ".\n" )
        std::getline( structfile, line );
        std::istringstream sstr( line );
        Crystal::Structure::t_Atom a;
        types::t_int type;
        sstr >> type;
        if( type != 1 or type != 2 ) continue;
        a.type = ( type == 1 ) ? -1.e0: 1.e0;
        sstr >> a.x[0] >> a.x[1]  >> a.x[2];
        a.freeze = Crystal::Structure::t_Atom::FREEZE_NONE;
        a.site = 0;
      }
      __ASSERT( structure.atoms.size() != N,
                "Could not read " << N << "atoms from file.\n" )
      
      // Adds structure to structure set.
      structures.push_back( structure );
      t_Configurations *confs = _symseps( structure );
      training.push_back( confs );
      delete confs; 
    }
    __CATCHCODE(, "Error while reading " << (path / _filename) << "\n" )
  }

}

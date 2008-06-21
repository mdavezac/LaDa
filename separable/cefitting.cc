//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string> 

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "cefitting.h" 

namespace Fitting
{
  void SepCeInterface::read( CE::SymSeparables &_symseps,
                             const std::string &_dir,
                             const std::string &_ldasdat )
  {
    try
    {
      boost::filesystem::path path( _dir );
      read_ldasdat( path, _ldasdat );
      std::vector< std::string > :: const_iterator i_name( names.begin() );
      std::vector< std::string > :: const_iterator i_name_end( names.end() );
      for(; i_name != i_name_end; ++i_name )
        read_structure( _symseps, path, *i_name );
    }
    __CATCHCODE(, "Error while reading training set.\n" )
  }


  void SepCeInterface :: read_ldasdat( const boost::filesystem::path &_path,
                                       const std::string &_ldasdat )
  {
    boost::filesystem::path fullpath = _path / _ldasdat;
    try
    {
      __ASSERT( boost::filesystem::exists( fullpath ), 
                "Cannot find " << fullpath << ".\n"   )
      std::ifstream ldas( fullpath.string().c_str(), std::ifstream::in );
      while( not ldas.eof() )
      {
        std::string name;
        types::t_real energy;
        ldas >> name >> energy;
        names.push_back( name );
        targets.push_back( energy );
        weight.push_back( 1 );
      }
    }
    __CATCHCODE(, "Error while reading " << fullpath << "\n" )
  }

  void SepCeInterface :: read_structure( CE::SymSeparables &_symseps,
                                         const boost::filesystem::path &_path, 
                                         const std::string &_filename )
  {
    boost::filesystem::path fullpath = _path / _filename;
    try
    {
      // Reads structure from structure file @ nrel.
      __ASSERT( boost::filesystem::exists( fullpath ),
                "Cannot find " << fullpath  << ".\n" )
      std::ifstream structfile( fullpath.string().c_str(), std::ifstream::in );
      std::string line;
      std::getline( structfile, line ); // name and inconsequential data.

      Crystal::Structure structure;
      types::t_int N;  // number of atoms;
      // cell 
      for(types::t_int i(0); i < 3; ++i )
      {
        __ASSERT( structfile.eof(),
                  "Reached unexpected end of file: " << fullpath << ".\n" )
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
        __ASSERT( structfile.eof(),
                  "Reached unexpected end of file: " << _filename << ".\n" )
        std::getline( structfile, line );
        std::istringstream sstr( line );
        Crystal::Structure::t_Atom a;
        types::t_int type;
        sstr >> type;
        if( type != 1 or type != 2 ) continue;
        a.type = ( type == 1 ) ? -1.e0: 1.e0;
        sstr >> a.pos.x[0] >> a.pos.x[1]  >> a.pos.x[2];
        a.freeze = Crystal::Structure::t_Atom::FREEZE_NONE;
        a.site = 0;
      }
      __ASSERT( structure.atoms.size() != N,
                "Could not read " << N << "atoms from " << fullpath << ".\n" )
      
      // Adds structure to structure set.
      structures.push_back( structure );
      t_Configurations *confs = _symseps.configurations( structure );
      training.push_back( *confs );
      delete confs; 
    }
    __CATCHCODE(, "Error while reading " << fullpath << "\n" )
  }

  std::pair<types::t_real, types::t_real> 
   SepCeInterface :: check( CE::Separables &_sep,
                            bool _verbose, bool _which ) const
    {
      types::t_real average(0);
      types::t_real maxerr(0);
      std::vector< t_Configurations > :: const_iterator i_train( training.begin() );
      std::vector< t_Configurations > :: const_iterator i_train_end( training.end() );
      std::vector< types::t_real > :: const_iterator i_weight( weight.begin() );
      std::vector< types::t_real > :: const_iterator i_target( targets.begin() );
      std::vector< std::string > :: const_iterator i_name( names.begin() );
      for( types::t_unsigned i=0; i_train != i_train_end ;
           ++i, ++i_train, ++i_weight, ++i_target, ++i_name )
      {
        // Checks whether fit or prediction.
        if( exclude.size() )
        {
          bool notfound
            = ( std::find( exclude.begin(), exclude.end(), i ) == exclude.end() );
          if( _which != notfound ) continue;
        }
        // Sums over all equivalent structure
        t_Configurations :: const_iterator i_conf( i_train->begin() );
        t_Configurations :: const_iterator i_conf_end( i_train->end() );
        types::t_real intermed(0);
        for(; i_conf != i_conf_end; ++i_conf )
          intermed +=  _sep( i_conf->first ) * ( i_conf->second );
        // Prints each structure if _verbose=true.
        if( _verbose )
          std::cout << "  structure: " << *i_name << "   "
                    << "Target: " << *i_target << " "
                    << "Separable: " << intermed << "   "
                    << "|Target-Separable| * weight: "
                    << std::abs( intermed - (*i_target) ) * (*i_weight) << "\n";
        intermed = std::abs( intermed - (*i_target) ) * (*i_weight);
        if( intermed > maxerr ) maxerr = intermed;
        average += intermed;
      }
      types::t_real div( training.size() );
      if( exclude.size() )
      {
        if( _which ) div  = (types::t_real ) exclude.size();
        else         div -= (types::t_real ) exclude.size();
      }
      average /= div;
      return std::pair< types::t_real, types::t_real>(average, maxerr);
    }

}
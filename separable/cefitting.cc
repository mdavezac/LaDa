//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string> 
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <opt/algorithms.h>

#include "cefitting.h" 

namespace Fitting
{
  void SepCeInterface::read( CE::SymSeparables &_symseps,
                             const std::string &_dir,
                             const std::string &_ldasdat,
                             bool _verbose )
  {
    try
    {
      boost::filesystem::path path( _dir );
      read_ldasdat( path, _ldasdat );
      {
        namespace bl = boost::lambda;
        std::for_each
        ( 
          names.begin(), names.end(),
          bl::bind( &SepCeInterface::read_structure, bl::var(*this),
                    bl::constant( _symseps ), bl::constant( path ), bl::_1 )
        );
        std::vector< Crystal::Structure > :: iterator i_s = structures.begin();
        std::vector< Crystal::Structure > :: iterator i_se = structures.end();
        std::vector< types::t_real > :: const_iterator i_t = targets.begin();
        for(; i_s != i_se; ++i_s, ++i_t ) i_s->energy = *i_t;
        if( not _verbose ) return;
        std::for_each
        ( 
          structures.begin(), structures.end(), 
          std::cout << bl::_1 << "\n" 
        );
      }
    }
    __CATCHCODE(, "Error while reading training set.\n" )
  }


  void SepCeInterface :: read_ldasdat( const boost::filesystem::path &_path,
                                       const std::string &_ldasdat )
  {
    boost::filesystem::path fullpath = _path / _ldasdat;
    try
    {
      __ASSERT( not boost::filesystem::exists( fullpath ), 
                "Cannot find " << fullpath << ".\n"   )
      std::ifstream ldas( fullpath.string().c_str(), std::ifstream::in );
      std::string line;
      while( std::getline( ldas, line ) )
      {
        const boost::regex re("^(\\s+)?(\\S+)\\s+(-?\\d+(\\.\\d+)?)");
        boost::match_results<std::string::const_iterator> what;
        if( not boost::regex_search( line, what, re ) ) continue;

        names.push_back( what.str(2) );
        targets.push_back( boost::lexical_cast<types::t_real>( what.str(3) ) );
        weight.push_back( 1 );
      }
    }
    __CATCHCODE(, "Error while reading " << fullpath << "\n" )
  }

  void SepCeInterface :: read_structure( const CE::SymSeparables &_symseps,
                                         const boost::filesystem::path &_path, 
                                         const std::string &_filename )
  {
    boost::filesystem::path fullpath = _path / _filename;
    try
    {
      // Reads structure from structure file @ nrel.
      __ASSERT( not boost::filesystem::exists( fullpath ),
                "Cannot find " << fullpath  << ".\n" )
      std::ifstream structfile( fullpath.string().c_str(), std::ifstream::in );
      std::string line;
      std::getline( structfile, line ); // name and inconsequential data.

      Crystal::Structure structure;
      structure.name = _filename;
      types::t_int N;  // number of atoms;
      std::getline( structfile, line ); // name and inconsequential data.
      std::istringstream sstr( line );
      sstr >> N;
      // cell 
      for(types::t_int i(0); i < 3; ++i )
      {
        __ASSERT( structfile.eof(),
                  "Reached unexpected end of file: " << fullpath << ".\n" )
        std::getline( structfile, line );
        sstr.str( line ); sstr.seekg (0, std::ios::beg); sstr.clear();
        sstr >> structure.cell.x[i][0]
             >> structure.cell.x[i][1]
             >> structure.cell.x[i][2];
      }
      structure.freeze = Crystal::Structure::FREEZE_NONE;
      // now atoms.
      types::t_int nfound(0);
      while( nfound < N and structfile.good() )
      {
        __ASSERT( structfile.eof(),
                  "Reached unexpected end of file: " << _filename << ".\n" )
        std::getline( structfile, line );
        sstr.str( line ); sstr.seekg (0, std::ios::beg); sstr.clear();
        Crystal::Structure::t_Atom a;
        types::t_int type;
        sstr >> type;
        ++nfound;
        if( type != 1 and type != 2 ) continue;
        a.type = ( type == 1 ) ? -1.e0: 1.e0;
        sstr >> a.pos.x[0] >> a.pos.x[1]  >> a.pos.x[2];
        a.freeze = Crystal::Structure::t_Atom::FREEZE_NONE;
        a.site = 0;
        structure.atoms.push_back(a);
      }
      __ASSERT( nfound != N,    "Could find only " << nfound << " of " 
                             << N << " atoms in " << fullpath << ".\n" )
      
      // Adds structure to structure set.
      structures.push_back( structure );
      training.resize( training.size() + 1 );
      t_Configurations &confs( training.back() );
      confs.clear();
      _symseps.configurations( structure, confs );
    }
    __CATCHCODE(, "Error while reading " << fullpath << "\n" )
  }

  std::pair<types::t_real, types::t_real> 
   SepCeInterface :: check( const CE::Separables &_sep,
                            bool _which, bool _verbose ) const
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
        {
          intermed +=  _sep.operator()( i_conf->first.begin(),
                                        i_conf->first.end()   ) * i_conf->second;
        }
        // Prints each structure if _verbose=true.
        if( _verbose )
          std::cout << "  structure: " << std::setw(30) << *i_name << "   "
                    << "Target: " << std::fixed << std::setw(8)
                           << std::setprecision(2) << *i_target + offset << " "
                    << "Separable: " << std::fixed << std::setw(8)
                           << std::setprecision(2) << intermed << "   "
                    << "|Target-Separable| * weight: "
                    << std::fixed << std::setw(10) << std::setprecision(3) 
                           << std::abs( intermed - (*i_target) - offset ) * (*i_weight) << "\n";
        intermed = std::abs( intermed - (*i_target) - offset ) * (*i_weight);
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

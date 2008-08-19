//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


#include <opt/random.h>
#include<boost/tokenizer.hpp>

#include "leave_many_out.h"


namespace Fitting
{

  void LeaveManyOut::extract_cmdl( boost::program_options::variables_map &_map )
  {
    try
    {
      namespace bl = boost::lambda;
      namespace po = boost::program_options;
      namespace fs = boost::filesystem;

      std::string line;
      boost::match_results<std::string::const_iterator> what;
      if( _map.count(cmdl_except) )
      {
        except.clear();
        line = _map[cmdl_except].as<std::string>();
        boost::char_separator<char> sep(", ");
        typedef boost::tokenizer< boost::char_separator<char> > t_Tokenizer;
        t_Tokenizer tok(line, sep);
        for( t_Tokenizer::iterator i_tok=tok.begin(); i_tok!=tok.end();++i_tok)
          except.push_back( boost::lexical_cast<types::t_unsigned>( *i_tok ) );
      }
      if( not ( _map.count(cmdl_set) or _map.count(cmdl_file) ) ) 
      {
        do_perform = false;
        sets.clear();
        return;
      }
      do_perform = true;

      if( _map.count(cmdl_file) ) 
        filename = _map[cmdl_file].as<std::string>();

      if( not _map.count(cmdl_set) ) 
      {
        read_sets();
        return;
      }

      line = _map[cmdl_set].as<std::string>();
      do_perform = true;
      
      const boost::regex re("(\\d+)\\s+(\\d+)");
      __DOASSERT( not boost::regex_search( line, what, re ),
                     "Could not make sense of --" << cmdl_set 
                  << "=\"" << line << "\".\n" )

      types::t_unsigned nbsets = boost::lexical_cast<types::t_unsigned>( what.str(1) );
      types::t_int perset = boost::lexical_cast<types::t_int>( what.str(2) );

      fs::path fullpath = filename;
      if( not fs::exists( fullpath ) ) goto out;
      
      read_sets();
      if( nb_sets ==  nbsets and perset == nb_perset ) return;
      
      sets.clear();
out:
      if( except.size() == 1 )
        std::cout << "   Structure " << except.front() << " is always included in fit.\n"; 
      else if( except.size() == 2 )
        std::cout << "   Structures " << except.front() << " and " 
                  << except.back() << " are always included in fit.\n"; 
      else if( except.size() >= 3 )
      {
        std::cout << "   Structures ";
        std::for_each( except.begin(), except.end()-1, std::cout << bl::_1 << ", " );
        std::cout << ", and " << except.back() << " are always included in fit.\n";
      }
      if( nb_sets ==  nbsets and perset == nb_perset ) return;
      nb_sets = nbsets;
      nb_perset = perset;
    }
    __CATCHCODE(, "Error while extracting LeaveManyOut command-line options.\n" );
  }

  void LeaveManyOut :: read_sets()
  {
    namespace fs = boost::filesystem;

    sets.clear();
    except.clear();

    nb_sets = 0;
    nb_perset = 0;
    fs::path fullpath = filename;
    __DOASSERT( not fs::exists( fullpath ),
                "Could not find file " << fullpath << ".\n" )

    std::ifstream ldas( fullpath.string().c_str(), std::ifstream::in );
    std::string line;
    std::cout << "Reading leave-many-out sets from: " << filename << ".\n";
    while( std::getline( ldas, line ) )
    {
      std::ostringstream sstr;
      std::vector< types::t_unsigned > set;
      boost::tokenizer<> tok(line);
      boost::tokenizer<>::iterator i_tok=tok.begin(); 

      if( (*i_tok)[0] == '#' ) continue;

      for(; i_tok!=tok.end(); ++i_tok )
      {
        set.push_back( boost::lexical_cast< types::t_unsigned >( *i_tok ) );
        sstr << set.back() << " ";
      }

      if( not set.size() ) continue;
      
      ++nb_sets;
      std::cout << "  _ " << sstr.str() << "\n";
      if( nb_perset == 0 ) nb_perset = set.size();
      else if( nb_perset != set.size() ) nb_perset = -1;
      sets.push_back(set);
    }
    __DOASSERT( sets.size() == 0, "Read 0 sets from input.\n" )
  }

  void LeaveManyOut :: create_sets( types::t_unsigned _tsize )
  {
    namespace fs = boost::filesystem;

    __DOASSERT( _tsize == 0, "Input error: number of input structures is 0.\n" )
    __DOASSERT( nb_sets == 0, "Input error: required creation of 0 sets.\n" )
    __DOASSERT( nb_perset == 0, "Input error: required creation of sets of size 0.\n" )

    fs::path fullpath = filename;
    std::ofstream setfile( fullpath.string().c_str(),
                           std::ofstream::out | std::ofstream::trunc );
    __ASSERT( not setfile.is_open(), "Could not open " << fullpath << " for writing.\n" )

    std::cout << "Creating leave-many-out sets ( saved to file: " 
              << fullpath << " ): \n";
    for(types::t_unsigned i(nb_sets); i > 0; --i )
    {
      std::ostringstream sstr;
      std::vector< types::t_unsigned > set;
      for(types::t_unsigned j(nb_perset); j > 0; --j )
      {
        types::t_unsigned r;
        do { r = opt::random::range(0, _tsize ); }
        while(     except.size() != 0
               and except.end() != std::find( except.begin(), except.end(), r ) );
        set.push_back( r );
        sstr << set.back() << " ";
      }
      std::cout << "  _ " << sstr.str() << "\n";
      setfile << sstr.str() << "\n";
      sets.push_back( set );
    }
  }

}

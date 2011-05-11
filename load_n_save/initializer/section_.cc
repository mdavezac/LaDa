//
//  Version: $Id: section_.cc 1227 2009-07-14 02:17:07Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "section_.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      void Section_::FlowControl :: did_parse(std::string const& _name )
      {
        if( not parsed_options )
          std::cerr << "Could not parse options in section " + _name + ".\n";
        if( not parsed_subsections )
          std::cerr << "Could not parse subsections in section " + _name + ".\n";
      
        found = parsed_options and parsed_subsections;
      }

      struct Fresh
      {
        bool operator()( tree::Section const& _sec ) const
          { return _sec.fresh; }
        bool operator()( tree::Option const& _sec ) const
          { return _sec.fresh; }
      };

      void Section_::print_unknowns_( std::string const& _name, tree::Section const& _sec ) const
      {
        typedef boost::filter_iterator<Fresh, tree::Section::const_iterator::subsection > t_itersec;
        t_itersec i_sec( _sec.subsections_begin(), _sec.subsections_end() );
        t_itersec const i_sec_end( Fresh(), _sec.subsections_end(), _sec.subsections_end() );
        if( i_sec != i_sec_end )
        {
          std::cerr << "Unknown sections in " + _name + ":";
          for(; i_sec != i_sec_end; ++i_sec )
            std::cerr << " " << i_sec->name;
          std::cerr << "\n";
        }

        typedef boost::filter_iterator<Fresh, tree::Section::const_iterator::option > t_iterop;
        t_iterop i_op( _sec.options_begin(), _sec.options_end() );
        t_iterop const i_op_end( _sec.options_end(), _sec.options_end() );
        if( i_op != i_op_end )
        {
          std::cerr << "Unknown options in " + _name + ":";
          for(; i_op != i_op_end; ++i_op )
            std::cerr << " " << i_op->name;
          std::cerr << "\n";
        }
      }
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa



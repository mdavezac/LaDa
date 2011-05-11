//
//  Version: $Id: printer.cc 1226 2009-07-13 06:28:01Z davezac $
//


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../tree/tree.h"

#include "printer.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xml
    {
      std::ostream& print( std::ostream& _stream, 
                           tree::Base const &_base,
                           std::string const &_indent, 
                           std::string const &_char )
      {
        tree::Base :: const_iterator :: subsection i_begin( _base.subsections_begin() );
        tree::Base :: const_iterator :: subsection i_end( _base.subsections_end() );
        for(; i_begin != i_end; ++i_begin ) print( _stream, *i_begin, _indent, _char );
        return _stream;
      }


      std::ostream& print( std::ostream& _stream, 
                           tree::Section const &_sec,
                           std::string const &_indent,
                           std::string const &_char )
      {
        std::string const str( _indent + "<" + _sec.name );
        size_t const s( str.size() + 1 );
        _stream << str;
        { // print options.
          tree::Section :: const_iterator :: option i_begin( _sec.options_begin() );
          tree::Section :: const_iterator :: option i_end( _sec.options_end() );
          for( size_t i(s); i_begin != i_end; ++i_begin )
          {
            std::string op( i_begin->name );
            if( i_begin->name.size() > 0 )
              op += "=\"" + i_begin->value + "\"";
            if( i + op.size()  > 100 and s < 100 )
            {
              i = s;
              _stream << "\n";
              for( size_t j(s); j > 0; --j )
                _stream << " ";
            }
            _stream << " " << op;
          }
        }
        tree::Section :: const_iterator :: subsection i_begin( _sec.subsections_begin() );
        tree::Section :: const_iterator :: subsection const i_end( _sec.subsections_end() );
        if( i_begin == i_end ) 
        {
          _stream << "/>\n"; 
          return _stream;
        }
        _stream << ">\n"; 
        for(; i_begin != i_end; ++i_begin ) print( _stream, *i_begin, _indent + _char, _char );
        _stream << _indent << "</" << _sec.name << ">\n";
        return _stream;
      }
    } // namespace xml
  } // namespace load_n_save
} // namespace LaDa


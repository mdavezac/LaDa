//
//  Version: $Id: binary.cc 1266 2009-08-10 05:01:26Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <iterator>

#include <opt/debug.h>
#include "binary.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace sequencer
    {
      bool is_all_ands( Binary const &_a ) 
      {
        size_t i,j;
        Binary::const_iterator i_first( _a.begin() );
        Binary::const_iterator const i_end( _a.end() );
        for(; i_first != i_end; ++i_first )
          if( *i_first == first_object or *i_first == second_object ) continue;
          else if( *i_first == start_group )
          {
            ++i_first; 
            LADA_DOASSERT( i_first != i_end, "Unexpected end of container.\n" )
            details::find_group_end_impl( i_first, i_end, i, j );
          }
          else if( *i_first == or_ ) return false;
#         ifdef _LADADEBUG
           else if( *i_first == end_group )
             __THROW_ERROR( "Unexpected end of group.\n" );
#         endif
        return true;
      }

      void Binary :: transform_( std::string const &_str )
      {
        try
        {
          std::string :: const_iterator i_first = _str.begin();
          std::string :: const_iterator const i_end = _str.end();
          for(; i_first != i_end; ++i_first )
            if( *i_first == '1' )      std::list<tags>::push_back(first_object);
            else if( *i_first == '2' ) std::list<tags>::push_back(second_object);
            else if( *i_first == '(' ) std::list<tags>::push_back(start_group);
            else if( *i_first == ')' ) std::list<tags>::push_back(end_group);
            else if( *i_first == '|' ) std::list<tags>::push_back(or_);
            else if( *i_first == '&' ) continue;
            else
            { 
#             ifdef _LADA_DEBUG 
                std::cerr <<  "Unknown " << *i_first
                          << " while constructing Sequence object from string.\n";
#             endif
              std::list<tags>::clear();
              break;
            }
        }
        catch(...)
        {
          std::list<tags> :: clear();
#         ifdef _LADA_DEBUG 
            std::cerr <<  "Caught throw while constructing Sequence object from string.\n";
#         endif
        }
      }

      void operator&=( Binary &_a, Binary const &_b )
      {
        size_t const na( _a.size() );
        size_t const nb( _b.size() );

        if( na == 0 and nb == 0 ) return;
        if( na == 0 ) { _a = _b; return; }
        if( nb == 0 ) return;

        bool const a_isall_ands( is_all_ands( _a ) );
        bool const b_isall_ands( is_all_ands( _b ) );

        if( not a_isall_ands )
        {
          _a.push_front(start_group);
          _a.push_back(end_group);
        }
        if( not b_isall_ands ) _a.push_back(start_group);
        std::copy( _b.begin(), _b.end(), std::back_inserter(_a) );
        if( not b_isall_ands ) _a.push_back(end_group);

      }

#     ifdef _LADADEBUG
        std::ostream& operator<<( std::ostream &_stream, Binary const& _s )
        {
          Binary::const_iterator i_first( _s.begin() );
          Binary::const_iterator const i_end( _s.end() );
          for(; i_first != i_end; ++i_first )
            switch(*i_first)
            {
              case first_object: _stream << '1'; break;
              case second_object: _stream << '2'; break;
              case start_group: _stream << '('; break;
              case end_group: _stream << ')'; break;
              case or_: _stream << '|'; break;
            }
          return _stream;
        }
#     endif

    } // namespace sequencer.
  } // namespace load_n_save

} // namespace LaDa



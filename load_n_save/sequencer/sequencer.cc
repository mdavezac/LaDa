//
//  Version: $Id: sequencer.cc 1266 2009-08-10 05:01:26Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <iterator>

#include <opt/debug.h>
#include "sequencer.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace sequencer
    {
      bool is_all_ands( Sequence const &_a ) 
      {
        size_t i;
        Sequence::const_iterator i_first( _a.begin() );
        Sequence::const_iterator const i_end( _a.end() );
        for(; i_first != i_end; ++i_first )
          if( *i_first == object ) continue;
          else if( *i_first == start_group )
          {
            ++i_first; 
            LADA_DOASSERT( i_first != i_end, "Unexpected end of container.\n" )
            details::find_group_end_impl( i_first, i_end, i );
          }
          else if( *i_first == or_ ) return false;
#         ifdef _LADADEBUG
           else if( *i_first == end_group )
             __THROW_ERROR( "Unexpected end of group.\n" );
#         endif
        return true;
      }

      void Sequence :: transform_( std::string const &_str )
      {
        try
        {
          std::string :: const_iterator i_first = _str.begin();
          std::string :: const_iterator const i_end = _str.end();
          for(; i_first != i_end; ++i_first )
            if( *i_first == 'o' )      std::list<tags>::push_back(object);
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

      void operator&=( Sequence &_a, Sequence const &_b )
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
        std::ostream& operator<<( std::ostream &_stream, Sequence const& _s )
        {
          Sequence::const_iterator i_first( _s.begin() );
          Sequence::const_iterator const i_end( _s.end() );
          for(; i_first != i_end; ++i_first )
            switch(*i_first)
            {
              case object: _stream << 'o'; break;
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



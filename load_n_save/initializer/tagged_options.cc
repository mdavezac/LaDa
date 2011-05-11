//
//  Version: $Id: tagged_options.cc 1226 2009-07-13 06:28:01Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tagged_options.h"

namespace LaDa 
{
  namespace load_n_save
  {

    namespace details
    {
      template< class T1 > 
        bool do_options( T1 &_first, T1 const &_end, tree::Section const & _sec ) 
        {
          bool found_all = true;
          ++_first;
          for(; _first != _end; ++_first )
          {
            if( *_first == "(" )  
            {
              found_all &= do_options( _first, _end, _sec ); 
              continue;
            }
            if( *_first == ") or (" )  // between two groups
            {
              if( not found_all ) continue;
              break;
            }
            if( *_first == ")" ) break; // end of groups.
            found_all &=  _sec.options_begin( *_first ) != _sec.options_end( *_first );
          } 
          return found_all;
        }
    } // namespace details

    bool TaggedOptions :: operator()( tree::Section const & _sec ) const
    {
      if( not ( (bool) options_ ) ) return true;
      TaggedOptions::t_Options :: const_iterator i_op = options_->begin();
      TaggedOptions::t_Options :: const_iterator const i_end = options_->end();
      if( i_op == i_end ) return true;
      for(; i_op != i_end; ++i_op )
      {
        if( *i_op == "(" )  // marks beginning of groups
          details::do_options( i_op, i_end, _sec );
        else if( _sec.options_begin( *i_op ) == _sec.options_end( *i_op ) ) return false;
      }
      return true;
    }
    

    TaggedOptions TaggedOptions :: unknown( tree::Section const & _sec ) const
    {
      TaggedOptions result(*this);
      result.options_.reset( new t_Options );
      result.tag_ = tags::option::default_;
      tree::Section::const_iterator::option i_op = _sec.options_begin();
      tree::Section::const_iterator::option const i_end = _sec.options_end();
      if( i_op == i_end ) return result;
      for(; i_op != i_end; ++i_op )
        if( not ( (bool) options_ ) ) 
          result.options_->push_back( i_op->name );
        else if( 
                 options_->end() == std::find
                                     ( 
                                       options_->begin(), 
                                       options_->end(), 
                                       i_op->name
                                     ) 
               )
          result.options_->push_back( i_op->name );
      return result;
    }

    std::ostream& operator<<( std::ostream& _stream, TaggedOptions const &_idops )
    {
      if( not ( (bool) _idops.options_ ) ) return _stream << "none.\n";
      if( _idops.options_->size() == 0 ) return _stream << "none.\n";
      TaggedOptions::t_Options :: const_iterator i_op = _idops.options_->begin();
      TaggedOptions::t_Options :: const_iterator const i_end = _idops.options_->end();
      bool is_first_of_op(true);
      for(; i_op != i_end; ++i_op )
      {
        bool with_comma(true);
        if( *i_op  == "(" or *i_op == ") or (" ) 
        {
          is_first_of_op = true;
          with_comma = false;
        }
        else if( *i_op == ")" )
          with_comma = false;
        else if ( is_first_of_op )
        {
          is_first_of_op = false;
          with_comma = false;
        }
        if( with_comma ) _stream << ", ";
        _stream << *i_op;
      }
      return _stream << ".\n";
    }
  } // namespace load_n_save

} // namespace LaDa



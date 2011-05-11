//
//  Version: $Id: tagged_section.cc 1201 2009-06-22 05:26:09Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tagged_sections.h"

namespace LaDa 
{
  namespace load_n_save
  {
    bool TaggedSections :: operator()( tree::Section const & _sec ) const
    {
      TaggedSections::t_Sections :: const_iterator i_op = sections_.begin();
      TaggedSections::t_Sections :: const_iterator const i_end = sections_.end();
      if( i_op == i_end ) return true;
      for(; i_op != i_end; ++i_op )
        if( _sec.sections_begin( *i_op ) == _sec.sections_end( *i_op ) ) return false;
      return true;
    }
    //! Dumps id options to a stream.
    std::ostream& operator<<( std::ostream& _stream, TaggedSections const &_idops )
    {
      if( _idops.sections_.size() == 0 ) return _stream << "none.\n";
      TaggedSections::t_Sections :: const_iterator i_op = _idops.sections_.begin();
      TaggedSections::t_Sections :: const_iterator const i_end = --_idops.sections_.end();
      for(; i_op != i_end; ++i_op ) _stream << *i_op << ", ";
      return _stream << *i_op << ".\n";
    }
  } // namespace load_n_save

} // namespace LaDa



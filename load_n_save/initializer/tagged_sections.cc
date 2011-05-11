//
//  Version: $Id: tagged_sections.cc 1226 2009-07-13 06:28:01Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include "tagged_sections.h"

namespace LaDa 
{
  namespace load_n_save
  {
    bool TaggedSections :: operator()( tree::Section const & _sec ) const
    {
      if( not ( (bool) sections_ ) ) return true;
      TaggedSections::t_Sections :: const_iterator i_op = sections_->begin();
      TaggedSections::t_Sections :: const_iterator const i_end = sections_->end();
      if( i_op == i_end ) return true;
      for(; i_op != i_end; ++i_op )
        if( _sec.subsections_begin( *i_op ) == _sec.subsections_end( *i_op ) ) return false;
      return true;
    }

    TaggedSections TaggedSections :: unknown( tree::Section const & _sec ) const
    {
      TaggedSections result(*this);
      result.sections_.reset( new t_Sections );
      result.tag_ = tags::section::default_;
      tree::Section::const_iterator::subsection i_op = _sec.subsections_begin();
      tree::Section::const_iterator::subsection const i_end = _sec.subsections_end();
      if( i_op == i_end ) return result;
      for(; i_op != i_end; ++i_op )
        if( not ( (bool) sections_ ) ) 
          result.sections_->push_back( i_op->name );
        else if( 
                 sections_->end() == std::find
                                     ( 
                                       sections_->begin(), 
                                       sections_->end(), 
                                       i_op->name 
                                      ) 
               )
          result.sections_->push_back( i_op->name );
      return result;
    }
    
    std::ostream& operator<<( std::ostream& _stream, TaggedSections const &_idops )
    {
      if( not ( (bool) _idops.sections_ ) ) return _stream << "none.\n";
      if( _idops.sections_->size() == 0 ) return _stream << "none.\n";
      TaggedSections::t_Sections :: const_iterator i_op = _idops.sections_->begin();
      TaggedSections::t_Sections :: const_iterator const i_end = --_idops.sections_->end();
      for(; i_op != i_end; ++i_op ) _stream << *i_op << ", ";
      return _stream << *i_op << ".\n";
    }
  } // namespace load_n_save

} // namespace LaDa



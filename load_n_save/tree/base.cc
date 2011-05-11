//
//  Version: $Id: base.cc 1227 2009-07-14 02:17:07Z davezac $
//


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "base.h"
#include "section.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace tree
    {
      Base::t_Subsections::iterator Base::subsections_begin() 
      { 
        if( not subsections_ ) 
        {
          static iterator::subsection a;
          return a;
        }
        return iterator::subsection( subsections_->begin(), subsections_->end() );
      }
      Base::t_Subsections::const_iterator Base::subsections_begin() const
      { 
        if( not subsections_ ) 
        {
          static const_iterator::subsection a;
          return a;
        }
        return const_iterator::subsection( subsections_->begin(), subsections_->end() );
      }
      
      void Base::insert_subsection( Section const &_value )
      {
        if( not subsections_ ) subsections_.reset( new t_Subsections );
        subsections_->insert(_value);
      }
      
      void Base::copy_to( Base& _c ) const
      {
        if( subsections_ ) _c.subsections_.reset( new t_Subsections( *subsections_ ) );
        else if( _c.subsections_ ) 
        {
          boost::shared_ptr<t_Subsections> s;
          _c.subsections_.swap(s);
        }
      }

    } // namespace parser
  } // namespace load_n_save
} // namespace LaDa


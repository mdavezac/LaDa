//
//  Version: $Id: section.h 1266 2009-08-10 05:01:26Z davezac $
//

#ifndef _LADA_LOADNSAVE_PARSER_BASE_SECTION_H_
#define _LADA_LOADNSAVE_PARSER_BASE_SECTION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../string_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      // Forward declaration.
      class Section;
      // Forward declaration.
      class regular_data;
    }
    namespace parser_base
    {
     
      //! Abstract base class.
      struct Section
      {
        //! Virtual destructor.
        virtual ~Section() {}
        //! Double dispatch.
        virtual bool regular( xpr::Section const &_sec, xpr::regular_data const& ) const = 0;
        //! Parses content.
        virtual bool content( xpr::Section const&, t_String const& _n = "" ) const = 0;
        //! Parses subsections.
        virtual bool operator&( xpr::Section const& ) const = 0;
      };

    } // namespace parser_base
  } // namespace load_n_save

} // namespace LaDa


#endif

//
//  Version: $Id: tags.h 1227 2009-07-14 02:17:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_TAGS_H_
#define _LADA_LOADNSAVE_GRAMMAR_TAGS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa 
{
  namespace load_n_save
  {

    namespace tags
    {
      //! Scoped tags for sections.
      namespace section
      {
        //! Tags for sections.
        enum section
        {
          //! by default, all sections are optional.
          required = 1,
        };
        size_t const default_ = 0;
      }

      //! Scoped tags for sections.
      namespace option
      {
        // Tags for sections.
        enum option
        {
          //! by default, all options are optional.
          required = 1, 
          //! by default, all options should be unique.
          many     = 2,
          //! A option can be used to identify in conjunction with a section.
          id       = 4,  
        };
        size_t const default_ = 0;
      }
    } // namespace Tags

  } // namespace load_n_save

} // namespace LaDa


#endif

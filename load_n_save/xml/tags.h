//
//  Version: $Id: tags.h 1167 2009-06-07 23:22:42Z davezac $
//

#ifndef _LADA_LNS_XML_TAGS_H_
#define _LADA_LNS_XML_TAGS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xml
    {
      //! Offers scoped enums.
      namespace tags
      {
        //! Reading Tags.
        enum read 
        {
          READ = 1,
          VALIDATE = 2,
          READ_N_VALIDATE = 3,
        };
      };
    } // namespace xml
  } // namespace load_n_save
} // namespace LaDa

#endif

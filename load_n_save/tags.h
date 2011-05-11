//
//  Version: $Id: tags.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LNS_TAGS_H_
#define _LADA_LNS_TAGS_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace LaDa
{
  namespace load_n_save
  {
    enum tags
    {
      required = 1, //!< Required by a section.
      id       = 2  //!< Helps identify a section. Options only.
    };
  };
} // namespace LaDa


#endif 

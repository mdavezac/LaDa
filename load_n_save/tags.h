#ifndef LADA_LNS_TAGS_H
#define LADA_LNS_TAGS_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace LaDa
{
  namespace load_n_save
  {
    enum tags
    {
      unavailable = -1, //!< Reserved.
      optional    =  0, //!< Optional section or option.
      required    =  1, //!< Required section or option.
      idoption    =  2  //!< Helps identify a section. Options only.
    };
  };
} // namespace LaDa


#endif 

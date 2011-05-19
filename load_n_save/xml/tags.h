#ifndef LADA_LNS_XML_TAGS_H
#define LADA_LNS_XML_TAGS_H

#include "LaDaConfig.h"

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

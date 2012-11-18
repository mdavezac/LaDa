#ifndef LADA_CRYSTAL_PRIMITIVE_H
#define LADA_CRYSTAL_PRIMITIVE_H
#include "LaDaConfig.h"

#include "structure/structure.h"


namespace LaDa
{
  namespace crystal 
  {
#   ifdef LADA_CRYSTAL_MODULE
      //! Returns the primitive unit structure. 
      static Structure primitive(Structure const &_structure, types::t_real _tolerance = -1e0);
      //! Returns True if the input is primitive.
      static bool is_primitive(Structure const &_structure, types::t_real _tolerance = -1e0);
#   else
      //! Returns the primitive unit structure. 
      inline Structure primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
        { return (*(Structure(*)(Structure const&, types::t_real))
                  api_capsule[12])(_structure, _tolerance); }
      //! Returns True if the input is primitive.
      inline bool is_primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
        { return (*(bool(*)(Structure const&, types::t_real))
                  api_capsule[13])(_structure, _tolerance); }
#   endif
  } // namespace crystal
} // namespace LaDa
#endif

#ifndef LADA_CRYSTAL_PRIMITIVE_H
#define LADA_CRYSTAL_PRIMITIVE_H
#include "LaDaConfig.h"

#include "structure/structure.h"


namespace LaDa
{
  namespace crystal 
  {
    namespace
    {
      //! Returns the primitive unit structure. 
      Structure primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
#     ifdef LADA_CRYSTAL_MODULE
          ;
#     else
          { return (*(Structure(*)(Structure const&, types::t_real))
                    api_capsule[12])(_structure, _tolerance); }
#     endif
      //! Returns True if the input is primitive.
      bool is_primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
#     ifdef LADA_CRYSTAL_MODULE
          ;
#     else
        { return (*(bool(*)(Structure const&, types::t_real))
                  api_capsule[13])(_structure, _tolerance); }
#     endif
    } // anonymous namespace
  } // namespace crystal
} // namespace LaDa
#endif

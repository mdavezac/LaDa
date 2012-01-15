#ifndef LADA_CRYSTAL_PRIMITIVE_H
#define LADA_CRYSTAL_PRIMITIVE_H
#include "LaDaConfig.h"

#include "structure.h"


namespace LaDa
{
  namespace crystal 
  {
    //! Returns the primitive unit structure. 
    Structure primitive(Structure const &_structure, types::t_real _tolerance = -1e0);
    //! Returns True if the input is primitive.
    bool is_primitive(Structure const &_structure, types::t_real _tolerance = -1e0);
  } // namespace crystal
} // namespace LaDa
#endif

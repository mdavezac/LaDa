#ifndef _LADA_PYTHON_VFF_HPP_
#define _LADA_PYTHON_VFF_HPP_

#include "LaDaConfig.h"

#include <crystal/structure.h>

namespace LaDa
{
  namespace python
  {
    void expose_vff();
    void expose_layeredvff();

    template<class T>
      void print_escan_input( T &_self, const std::string &_path,
                              Crystal::TStructure<std::string> const &_str) 
      {
        _self.set_structure(_str);
        _self.init(true, false);
        _self.print_escan_input( _path ); 
      }
  }
} // namespace LaDa
#endif

//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_VFF_HPP_
#define _LADA_PYTHON_VFF_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <utility>

#include "../functional.h"
#include "../layered.h"
#include "../va.h"

namespace LaDa
{
  namespace python
  {
    void expose_vff();
    void expose_layeredvff();

    //! Vff "functional" for python.
    typedef std::pair<Vff::VABase<Vff::Functional>, Crystal::Structure>  t_Vff;
    //! Vff "functional" for python.
    typedef std::pair<Vff::VABase<Vff::Layered>, Crystal::Structure>  t_LayeredVff;

    template<class T>
      void print_escan_input( T &_self, const std::string &_path,
                              Crystal::TStructure<std::string> const &_str) 
      {
        Crystal::convert_string_to_real_structure(_str, _self.second);
        _self.first.init(true, false);
        _self.first.print_escan_input( _path ); 
      }
  }
} // namespace LaDa
#endif

//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_ESCAN_HPP_
#define _LADA_PYTHON_ESCAN_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/structure.h>
#include "../interface.h"

namespace LaDa
{
  namespace python
  {
    void expose_escan();
    void expose_genpot();
    bool create_directory(Pescan::Interface::t_Path const &_path);
    void set_scale(Pescan::Interface &_interface, Crystal::TStructure<std::string> const &_str);
  }
} // namespace LaDa
#endif

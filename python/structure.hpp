//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/structure.h>

namespace PythonLaDa
{
  namespace XML
  {
    template<> std::string nodename<Crystal::Structure>();
  }
  void export_structure();
}

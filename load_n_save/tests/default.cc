#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/fusion/tuple/tuple_tie.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../save/save.h"
#include "../action/enum.h"
#include "../exceptions.h"
#include "../action/fusion.h"

#include "math/serialize.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    float f;
    A() : f(0e0) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          namespace lns = LaDa :: load_n_save;
          return _ar & lns::section("A") << lns::option("f", lns::action=f, lns::default_=3.1416);
        }
  };

}

int main(int argc, char *argv[]) 
{
  here::A a;
  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > xml;
  lns::load::Load loader;
  bool result;

  a.f = 0e0;
  xml = lns::xml::parse("<A/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(abs(a.f-3.1416)<1e-12, "default did not work.");

  a.f = 0e0;
  xml = lns::xml::parse("<A f=\"-2.71828\" />");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(abs(a.f+2.71828)<1e-12, "default did not work.");

  return 0;
}

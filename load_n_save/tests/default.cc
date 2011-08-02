#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/fusion/tuple/tuple_tie.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../exceptions.h"

#include "math/serialize.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

#if LADA_TEST_TYPE == 0
#  define LADA_TYPE float
#  define LADA_INIT 0e0
#  define LADA_DEFAULT 3.1416
#  define LADA_LOAD std::abs(a.f-3.1416) < 1e-6
#  define LADA_XML "<A f=\"2.71828\"/>"
#  define LADA_RELOAD std::abs(a.f-2.71828) < 1e-6
#elif LADA_TEST_TYPE == 1
#  define LADA_TYPE std::string
#  define LADA_INIT ""
#  define LADA_DEFAULT "hello"
#  define LADA_LOAD a.f == "hello"
#  define LADA_XML "<A f=\"goodbye\"/>"
#  define LADA_RELOAD a.f == "goodbye"
#elif LADA_TEST_TYPE == 2
#  define LADA_TYPE math::rVector3d
#  define LADA_INIT math::rVector3d(0,0,0)
#  define LADA_DEFAULT math::rVector3d(0.5,0.5,0.5)
#  define LADA_LOAD (a.f - LADA_DEFAULT).squaredNorm() < 1e-12
#  define LADA_XML "<A f=\"0.2 -0.1 0.2\"/>"
#  define LADA_RELOAD (a.f - math::rVector3d(0.2, -0.1, 0.2)).squaredNorm() < 1e-12
#endif

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    LADA_TYPE f;
    A() : f(LADA_INIT) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          namespace lns = LaDa :: load_n_save;
          return _ar & lns::section("A") << lns::option("f", lns::action=f, lns::default_=LADA_DEFAULT);
        }
  };

}

int main(int argc, char *argv[]) 
{
  here::A a;
  boost::shared_ptr< lns::tree::Base > xml;
  lns::load::Load loader;
  bool result;

  a.f = LADA_INIT;
  xml = lns::xml::parse("<A/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.\n");
  LADA_DOASSERT(LADA_LOAD, "default did not work.\n");
  a.f = LADA_INIT;
  xml = lns::xml::parse(LADA_XML);
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.\n");
  LADA_DOASSERT(LADA_RELOAD, "default did not work.\n");
  return 0;
}

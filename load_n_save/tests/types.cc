#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"

#include "../save/save.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

#if LADA_INC_TYPE == 0
#  define LADA_TYPE int
#  define LADA_INIT 1
#  define LADA_XML "<A attr=\"1\"/>"
#  define LADA_TWICE a.type == 2
#  define LADA_RELOAD a.type == 1
#elif LADA_INC_TYPE == 1
#  define LADA_TYPE std::string
#  define LADA_INIT "hello"
#  define LADA_XML "<A attr=\"hello\"/>"
#  define LADA_TWICE a.type == "hellohello"
#  define LADA_RELOAD a.type == "hello"
#elif LADA_INC_TYPE == 2
#  include <math/serialize.h>
#  define LADA_TYPE math::rVector3d
#  define LADA_INIT math::rVector3d(1.5,0.5,1.5)
#  define LADA_XML "<A attr=\"1.5 0.5 1.5\"/>"
#  define LADA_TWICE (a.type-math::rVector3d(3,1,3)).squaredNorm() < 1e-12
#  define LADA_RELOAD (a.type-math::rVector3d(1.5,0.5,1.5)).squaredNorm() < 1e-12
#elif LADA_INC_TYPE == 3
#  include <math/serialize.h>
#  define LADA_TYPE LaDa::types::t_real
#  define LADA_INIT -3.1416
#  define LADA_XML "<A attr=\"-3.1416\"/>"
#  define LADA_TWICE abs(a.type+6.2832) < 1e-12
#  define LADA_RELOAD (a.type+3.1416) < 1e-12
#endif

namespace here
{
  class A
  {
#   ifdef LADA_INTERNAL
      friend class load_n_save::access;
#   endif
      
    public:
    LADA_TYPE type;
    A() { type = LADA_INIT; }
    LADA_TYPE get() const { return type;}
    void twice() { type += type; }
    private:
#   ifdef LADA_INTERNAL
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar & ( lns::section("A") << lns::option("attr", lns::action=type) );
        }
#   endif
  };
# ifndef LADA_INTERNAL
  template<class T_ARCHIVE>
    bool lns_access(T_ARCHIVE &_ar, A &_a, lns::version_type const _version)
    {
      return _ar & ( lns::section("A") << lns::option("attr", lns::action=_a.type) );
    }
# endif
}
int main(int argc, char *argv[]) 
{
  here::A a;
  lns::save::Save saver;
  // print to xml
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML, "Error in output.")

  // reload from xml
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( sstr.str() ) );
  a.twice();
  LADA_DOASSERT(LADA_TWICE, "error in twice.")
  lns::load::Load loader;
  bool result = loader( *xml, lns::ext(a) );
  LADA_DOASSERT(LADA_RELOAD, "error when reloading")

  return 0;
}

#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../action/set.h"
#include "../save/save.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

#if LADA_TEST_TYPE == 0
#  define LADA_TYPE int
#  define LADA_INIT type.insert(0); type.insert(5); type.insert(-6)
#  define LADA_XML "<A attr=\"-6 0 5\"/>"
#elif LADA_TEST_TYPE == 1
#  define LADA_TYPE std::string
#  define LADA_INIT type.insert("H"); type.insert("He"); type.insert("Li")
#  define LADA_XML "<A attr=\"H He Li\"/>"
#elif LADA_TEST_TYPE == 2
#  include <math/serialize.h>
#  define LADA_TYPE types::t_real
#  define LADA_INIT type.insert(0); type.insert(0.5); type.insert(-1.5)
#  define LADA_XML "<A attr=\"-1.5 0 0.5\"/>"
#endif

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    std::set<LADA_TYPE> type;
    A() { LADA_INIT; }
    private:
    template<class T_ARCHIVE>
      bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
      {
        return _ar & ( lns::section("A") << lns::option("attr", lns::action=type) );
      }
  };
}
int main(int argc, char *argv[]) 
{
  here::A a, b;
  lns::save::Save saver;
  // print to xml
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML, "Error in output:\n" + xmlstring + "\n")

  // reload from xml
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( sstr.str() ) );
  lns::load::Load loader;
  a.type.clear();
  LADA_DOASSERT(a.type != b.type, "error when clearing.\n")
  bool result = loader( *xml, lns::ext(a) );
  LADA_DOASSERT(a.type == b.type, "error when reloading.\n")

  return 0;
}

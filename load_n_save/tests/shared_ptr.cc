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

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    std::string type;
    A(std::string const &_hello = "hello") : type(_hello) {}
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
  boost::shared_ptr<here::A> a( new here::A() );
  lns::save::Save saver;
  // print to xml
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == "<A attr=\"hello\"/>", "Error in output.\n");

  // reload from xml
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse("<A attr=\"goodbye\"/>") );
  lns::load::Load loader;
  bool result = loader( *xml, lns::ext(a) );
  LADA_DOASSERT(a->type == "goodbye", "error when reloading.\n")

  boost::shared_ptr<here::A> null;
  result = true;
  xml = lns::xml::parse("<A attr=\"goodbye\"/>");
  try { loader.nocatch(*xml, lns::ext(null)); }
  catch(lns::error::empty_shared_ptr &_e) { result = false; }
  LADA_DOASSERT(result == false, "Did not fail as expected.\n");
  return 0;
}

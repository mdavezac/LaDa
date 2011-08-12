#include "LaDaConfig.h"

#define __PROGNAME__ "loadnsave"

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../save/save.h"
#include "../exceptions.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    int i;
    std::string hello;
    A() : i(1), hello("hello") {};
    void twice() { i += i; hello += hello; }
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar &
            (
              lns::section("A") << lns::option("i", lns::action=i, lns::tag=lns::optional) 
                                << lns::option("hello", lns::action=hello) 
            );
        }
  };
}

int main(int argc, char *argv[]) 
{
  std::string str1= "<A i=5 hello=\"goodbye\"/>";
  std::string str2= "<A hello=\"goodbye\"/>";
  std::string str3= "<A i=5/>";
  here::A a, b, c;
  lns::load::Load loader;
  bool result = true;
  boost::shared_ptr< lns::tree::Base > xml;


  xml = lns::xml::parse( str1 );
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "did not load")
  LADA_DOASSERT(a.i == 5, "Did not load i.")
  LADA_DOASSERT(a.hello == "goodbye", "Did not load hello.")

  xml = lns::xml::parse( str2 );
  LADA_DOASSERT(loader( *xml, lns::ext(b) ), "did not load")
  LADA_DOASSERT(b.i == 1, "Did not load i.")
  LADA_DOASSERT(b.hello == "goodbye", "Did not load hello.")

  xml = lns::xml::parse( str3 );
  try {  loader.nocatch( *xml, lns::ext(c) ); }
  catch(error::required_option_not_found &_e)  { result = false; }
  LADA_DOASSERT(result == false, "This should not have loaded.")

  return 0;
}

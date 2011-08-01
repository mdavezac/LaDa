#include "LaDaConfig.h"

#define __PROGNAME__ "loadnsave"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xpr/merge.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../save/save.h"
#include "../exceptions.h"
#include "../save/save.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    int i;
    A() : i(0) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar &
            (
              lns::section("A") << lns::option("i", lns::action=i)
            );
        }
  };

  class B 
  {
    friend class load_n_save::access;
    public:
    std::string hello;
    B() : hello("hello") {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar & ( lns::section("B") << lns::option("hello", lns::action=hello) );
        } 
  };

  class C : public A, public B
  {
    friend class load_n_save::access;
    public:
    float f;
    C() : A(), B(), f(-1.5) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar & lns::merge( lns::merge(static_cast<A*>(this), static_cast<B*>(this)),
                                   lns::section("C") << lns::option("f", lns::action=f) );
        }
  };
}

int main(int argc, char *argv[]) 
{
  here::C a;
  lns::load::Load loader;
  bool result = true;
  boost::shared_ptr< lns::tree::Base > xml;
  lns::save::Save saver;

  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == "<A i=\"0\" hello=\"hello\" f=\"-1.5\"/>", "Error in output.")
  LADA_DOASSERT(a.i == 0, "incorrect i.")
  LADA_DOASSERT(a.hello == "hello", "incorrect hello.")
  LADA_DOASSERT(std::abs(a.f+1.5)<1e-12, "incorrect f.")


  xml = lns::xml::parse("<A i=1 hello=\"goodbye\" f=\"0.5\"/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "did not load")
  LADA_DOASSERT(a.i == 1, "Did not load i.")
  LADA_DOASSERT(a.hello == "goodbye", "Did not load hello.")
  LADA_DOASSERT(std::abs(a.f-0.5)<1e-12, "Did not load f.")

  return 0;
}

#include "LaDaConfig.h"

#define __PROGNAME__ "loadnsave"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xpr/push_back.h"
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
    A(int _i = 0) : i(_i) {};
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
    std::vector<A> vec;
    int required;
    B(int _i=0) : required(_i) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          return _ar & lns::section("B") << lns::push_back(vec, required);
        }
  };
}

int main(int argc, char *argv[]) 
{
  here::B a;
  lns::load::Load loader;
  bool result = true;
  boost::shared_ptr< lns::tree::Base > xml;
  lns::save::Save saver;

  a.vec.push_back(here::A(0));
  a.vec.push_back(here::A(5));
  a.vec.push_back(here::A(-6));
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == "<B>  <A i=\"0\"/>  <A i=\"5\"/>  <A i=\"-6\"/></B>", "Error in output.")

  a.vec.clear();
  xml = lns::xml::parse(xmlstring);
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "did not load")
  LADA_DOASSERT(a.vec.size() == 3, "Did not load vec.")
  LADA_DOASSERT(a.vec[0].i == 0, "Did not load vec.")
  LADA_DOASSERT(a.vec[1].i == 5, "Did not load vec.")
  LADA_DOASSERT(a.vec[2].i == -6, "Did not load vec.")
  
  a.vec.clear();
  a.required = 2;
  xml = lns::xml::parse(xmlstring);
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "did not load")
  LADA_DOASSERT(a.vec.size() == 3, "Did not load vec.")
  LADA_DOASSERT(a.vec[0].i == 0, "Did not load vec.")
  LADA_DOASSERT(a.vec[1].i == 5, "Did not load vec.")
  LADA_DOASSERT(a.vec[2].i == -6, "Did not load vec.")

  a.vec.clear();
  a.required = 3;
  xml = lns::xml::parse(xmlstring);
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "did not load")
  LADA_DOASSERT(a.vec.size() == 3, "Did not load vec.")
  LADA_DOASSERT(a.vec[0].i == 0, "Did not load vec.")
  LADA_DOASSERT(a.vec[1].i == 5, "Did not load vec.")
  LADA_DOASSERT(a.vec[2].i == -6, "Did not load vec.")
 
  a.vec.clear();
  a.required = 4;
  result = true;
  xml = lns::xml::parse(xmlstring);
  try {loader.nocatch( *xml, lns::ext(a) ); }
  catch(error::too_few_sections &_e) { result = false; }
  LADA_DOASSERT(result == false, "Did not fail as expected.")

  return 0;
}

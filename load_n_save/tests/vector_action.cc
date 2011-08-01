#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../action/vector.h"
#include "../save/save.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

#if LADA_TEST_TYPE == 0
#  define LADA_TYPE int
#  define LADA_INIT type.push_back(0); type.push_back(5); type.push_back(-6)
#  define LADA_XML "<A attr=\"0 5 -6\"/>"
#  define LADA_TWICE  a.type[0] == 0 and a.type[1] == 10 and a.type[2] == -12
#  define LADA_RELOAD a.type[0] == 0 and a.type[1] == 5 and a.type[2] == -6
#elif LADA_TEST_TYPE == 1
#  define LADA_TYPE std::string
#  define LADA_INIT type.push_back("H"); type.push_back("He"); type.push_back("Li")
#  define LADA_XML "<A attr=\"H He Li\"/>"
#  define LADA_TWICE a.type[0] == "HH" and a.type[1] == "HeHe" and a.type[2] == "LiLi"
#  define LADA_RELOAD a.type[0] == "H" and a.type[1] == "He" and a.type[2] == "Li"
#elif LADA_TEST_TYPE == 2
#  include <math/serialize.h>
#  define LADA_TYPE types::t_real
#  define LADA_INIT type.push_back(0); type.push_back(0.5); type.push_back(-1.5)
#  define LADA_XML "<A attr=\"0 0.5 -1.5\"/>"
#  define LADA_TWICE     abs(a.type[0]) < 1e-12 \
                     and abs(a.type[1]-1.0) < 1e-12  \
                     and abs(a.type[2]+3.0) < 1e-12  
#  define LADA_RELOAD     abs(a.type[0]) < 1e-12 \
                      and abs(a.type[1]-0.5) < 1e-12  \
                      and abs(a.type[2]+1.5) < 1e-12  
#endif

namespace here
{
  class A
  {
    friend class load_n_save::access;
      
    public:
    std::vector<LADA_TYPE> type;
    A() { LADA_INIT; }
    void twice()
    { 
      std::vector<LADA_TYPE>::iterator i_first = type.begin();
      std::vector<LADA_TYPE>::iterator i_end = type.end();
      for(; i_first != i_end; ++i_first) *i_first += *i_first;
    }
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

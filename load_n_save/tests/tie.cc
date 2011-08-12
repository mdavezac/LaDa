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
    LaDa::math::rMatrix3d cell;
    A() 
    {
      for(size_t j(0); j < 3; ++j)
        for(size_t k(0); k < 3; ++k)
          cell(j,k) = j == k ? -0.5: 0.5;
    }
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          namespace bf = boost::fusion;
          namespace lns = LaDa :: load_n_save;
#         ifdef LADA_TIE
#            error LADA_TIE already defined.
#         endif
#         define LADA_TIE(i) bf::tie(cell(i,0), cell(i,1), cell(i,2))
#         ifdef LADA_TOE
#            error LADA_TOE already defined.
#         endif
#         define LADA_TOE(i) bf::tie(cell(0,i), cell(1,i), cell(2,i))
          lns::xpr::Section const seccell = lns::section("Cell") 
            << (
                    (    lns::option("r0", lns::tag=lns::required, lns::action=LADA_TIE(0))
                      && lns::option("r1", lns::tag=lns::required, lns::action=LADA_TIE(1))
                      && lns::option("r2", lns::tag=lns::required, lns::action=LADA_TIE(2)) )
                 || (    lns::option("a0", lns::tag=lns::required, lns::action=LADA_TOE(0))
                      && lns::option("a1", lns::tag=lns::required, lns::action=LADA_TOE(1))
                      && lns::option("a2", lns::tag=lns::required, lns::action=LADA_TOE(2))  )
               );
#         undef LADA_TIE
#         undef LADA_TOE
          return _ar & seccell;
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

  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring
      == "<Cell r0=\"-0.5 0.5 0.5\" r1=\"0.5 -0.5 0.5\" r2=\"0.5 0.5 -0.5\"/>", "Error in output.");
 
  a.cell = math::rMatrix3d::Zero();
  xml = lns::xml::parse("<Cell r0=\"-0.5 0.5 0.5\" r1=\"1.5 -0.5 0.5\" r2=\"0.5 0.5 -0.5\"/>" );
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  std::cout << a.cell << "\n";
  for(size_t i(0); i < 3; ++i)
    for(size_t j(0); j < 3; ++j)
      if(i == 1 and j == 0)
      { LADA_DOASSERT(abs(a.cell(i, j) - 1.5) < 1e-12, "Could not reload cell."); }
      else
      { LADA_DOASSERT(abs(a.cell(i, j) - (i==j?-0.5:0.5)) < 1e-12, "Could not reload cell."); }

  a.cell = math::rMatrix3d::Zero();
  xml = lns::xml::parse("<Cell a0=\"-0.5 0.5 0.5\" a1=\"1.5 -0.5 0.5\" a2=\"0.5 0.5 -0.5\"/>" );
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  for(size_t i(0); i < 3; ++i)
    for(size_t j(0); j < 3; ++j)
      if(j == 1 and i == 0)
      { LADA_DOASSERT(abs(a.cell(i, j) - 1.5) < 1e-12, "Could not reload cell."); }
      else
      { LADA_DOASSERT(abs(a.cell(i, j) - (i==j?-0.5:0.5)) < 1e-12, "Could not reload cell."); }

  result = true;
  a.cell = math::rMatrix3d::Zero();
  xml = lns::xml::parse("<Cell r0=\"-0.5 0.5 0.5\" a1=\"1.5 -0.5 0.5\" a2=\"0.5 0.5 -0.5\"/>" );
  try { loader.nocatch( *xml, lns::ext(a) ); }
  catch(error::required_option_not_found &_e) { result = false; }

  return 0;
}

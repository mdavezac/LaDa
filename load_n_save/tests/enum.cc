#include "LaDaConfig.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "../load/load.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../save/save.h"
#include "../action/enum.h"
#include "../exceptions.h"

using namespace LaDa;
namespace lns = LaDa :: load_n_save;

namespace here
{
  struct frozen
  {
    //! Tags to freeze cell coordinates.
    enum type
    {
      NONE       =  0,  //!< Freeze no atomic coordinate.
      X          =  1,  //!< Freeze x coordinate.
      Y          =  2,  //!< Freeze y coordinate.
      Z          =  4,  //!< Freeze z coordinate.
      T          =  8,  //!< Freeze type.
      CARTESIANS =  7,  //!< Freeze cartesians coordinates.
      ALL        =  15, //!< Freeze all.
    };
  };

  class A
  {
    friend class load_n_save::access;
      
    public:
    LaDa::types::t_unsigned i;
    bool inclusive;
    A(bool _inc = true) : i(frozen::X|frozen::Y), inclusive(_inc) {};
    private:
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version)
        {
          namespace lns = LaDa :: load_n_save;
          std::map<std::string, LaDa::types::t_unsigned> freeze_map;
          freeze_map["none"]      = frozen::NONE;
          freeze_map["x"]         = frozen::X;
          freeze_map["y"]         = frozen::Y;
          freeze_map["z"]         = frozen::Z;
          freeze_map["t"]         = frozen::T;
          freeze_map["cartesian"] = frozen::CARTESIANS;
          freeze_map["all"]       = frozen::ALL;
          return _ar & lns::section("A") 
                          << lns::option( "i",
                                          lns::action=lns::enum_(i, freeze_map, inclusive),
                                          lns::default_=frozen::NONE );
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

  a.i = here::frozen::X|here::frozen::Y|here::frozen::Z;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(a)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == "<A i=\"cartesian\"/>", "Error in output.");

  a.i = here::frozen::X|here::frozen::Y;
  other = saver(lns::ext(a));
  std::ostringstream sstr2;
  lns::xml::print(sstr2, *other);
  xmlstring = sstr2.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == "<A i=\"xy\"/>", "Error in output.");

  a.i = here::frozen::NONE;
  xml = lns::xml::parse( sstr.str() );
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(a.i == here::frozen::X|here::frozen::Y, "error when reloading");

  a.i = here::frozen::NONE;
  xml = lns::xml::parse("<A i=\"cartesian\"/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(a.i == here::frozen::X|here::frozen::Y|here::frozen::Z, "error when reloading");

  a.i = here::frozen::NONE;
  xml = lns::xml::parse("<A i=\"xz\"/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(a.i == here::frozen::X|here::frozen::Z, "error when reloading");

  result = true;
  a.i = here::frozen::NONE;
  xml = lns::xml::parse("<A i=\"cartsian\"/>");
  try { loader.nocatch( *xml, lns::ext(a) ); }
  catch(lns::error::option_parse_error &_e) { result = false; }
  LADA_DOASSERT(result == false, "did not throw as expected.");

  a.i = here::frozen::NONE;
  a.inclusive = false;
  xml = lns::xml::parse("<A i=\"xy\"/>");
  try { loader.nocatch( *xml, lns::ext(a) ); }
  catch(lns::error::option_parse_error &_e) { result = false; }
  LADA_DOASSERT(result == false, "did not throw as expected.");

  a.i = here::frozen::NONE;
  xml = lns::xml::parse("<A i=\"z\"/>");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(a.i == here::frozen::Z, "error when reloading");

  a.i = here::frozen::X;
  xml = lns::xml::parse("<A />");
  LADA_DOASSERT(loader( *xml, lns::ext(a) ), "Could not reload.");
  LADA_DOASSERT(a.i == here::frozen::NONE, "error when reloading");

  return 0;
}

#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include <load_n_save/load/load.h>
#include <load_n_save/xml/parser.h>
#include <load_n_save/xml/printer.h>
#include <load_n_save/save/save.h>
#include <opt/debug.h>
#include "../atom.h"

#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE std::string
#  define LADA_INIT_TYPE atom.type = "hello";
#  define LADA_CLEAR_TYPE atom.type = "";
#  define LADA_XML "<Atom pos=\"0.5 -0.5 0.5\" type=\"hello\" freeze=\"x\" site=\"-1\"/>"
#  define LADA_DOASSERT_TYPE LADA_DOASSERT(atom.type == "hello", "Wrong first type.\n");
#  define LADA_DOASSERT_SITE LADA_DOASSERT(atom.site == -1, "Wrong site index.\n");
#  define LADA_XML2 "<Atom pos=\"0.5 -0.5 0.5\" type=\"hello\" freeze=\"x\" site=\"1\"/>"
#  define LADA_DOASSERT_SITE2 LADA_DOASSERT(atom.site == 1, "Wrong site index.\n");
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE std::vector<std::string> 
#  define LADA_INIT_TYPE atom.type.push_back("Au"); atom.type.push_back("Pd");
#  define LADA_CLEAR_TYPE atom.type.clear();
#  define LADA_XML "<Atom pos=\"0.5 -0.5 0.5\" type=\"Au Pd\" freeze=\"x\" site=\"-1\"/>"
#  define LADA_DOASSERT_TYPE \
     LADA_DOASSERT(atom.type.size() == 2, "Wrong number of types.\n");\
     LADA_DOASSERT(atom.type[0] == "Au", "Wrong first type.\n");\
     LADA_DOASSERT(atom.type[1] == "Pd", "Wrong first type.\n");
#  define LADA_DOASSERT_SITE LADA_DOASSERT(atom.site == -1, "Wrong site index.\n");
#  define LADA_XML2 "<Atom pos=\"0.5 -0.5 0.5\" type=\"Au Pd\" freeze=\"x\" site=\"1\"/>"
#  define LADA_DOASSERT_SITE2 LADA_DOASSERT(atom.site == 1, "Wrong site index.\n");
#elif LADA_TEST_INCTYPE == 2
#  define LADA_TYPE std::set<std::string> 
#  define LADA_INIT_TYPE atom.type.insert("Au"); atom.type.insert("Pd");
#  define LADA_CLEAR_TYPE atom.type.clear();
#  define LADA_XML "<Atom pos=\"0.5 -0.5 0.5\" type=\"Au Pd\" freeze=\"x\" site=\"-1\"/>"
#  define LADA_DOASSERT_TYPE \
    { LADA_TYPE cmp; cmp.insert("Au"); cmp.insert("Pd");                \
      LADA_DOASSERT(atom.type.size() == 2, "Wrong number of types.\n");\
      LADA_DOASSERT(atom.type == cmp, "Wrong first type.\n"); }
#  define LADA_DOASSERT_SITE LADA_DOASSERT(atom.site == -1, "Wrong site index.\n");
#  define LADA_XML2 "<Atom pos=\"0.5 -0.5 0.5\" type=\"Au Pd\" freeze=\"x\" site=\"1\"/>"
#  define LADA_DOASSERT_SITE2 LADA_DOASSERT(atom.site == 1, "Wrong site index.\n");
#endif
using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  namespace lns = LaDa::load_n_save;
  Atom< LADA_TYPE > atom;
  atom.pos = math::rVector3d(0.5,-0.5,0.5);
  atom.freeze = Atom< LADA_TYPE >::frozen::X;
  LADA_INIT_TYPE;

  // print to xml
  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(atom)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML, "Error in output.");

  atom.pos = math::rVector3d(0,0,0);
  atom.freeze = Atom< LADA_TYPE >::frozen::NONE;
  LADA_CLEAR_TYPE;
  other = lns::xml::parse( sstr.str() );
  lns::load::Load loader;
  bool result = loader( *other, lns::ext(atom) );
  LADA_DOASSERT((atom.pos - math::rVector3d(0.5,-0.5,0.5)).squaredNorm() < 1e-12, "Could not reload pos.\n");
  LADA_DOASSERT(atom.freeze == Atom< LADA_TYPE >::frozen::X, "Could not reload freeze.\n");
  LADA_DOASSERT_TYPE
  LADA_DOASSERT_SITE

  // check site index stuff.
  atom.site = 1;
  other = saver(lns::ext(atom));
  std::ostringstream sstr2;
  lns::xml::print(sstr2, *other);
  xmlstring = sstr2.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML2, "Error in output.");
  
  atom.site = -1;
  other = lns::xml::parse( sstr2.str() );
  result = loader( *other, lns::ext(atom) );
  LADA_DOASSERT((atom.pos - math::rVector3d(0.5,-0.5,0.5)).squaredNorm() < 1e-12, "Could not reload pos.\n");
  LADA_DOASSERT(atom.freeze == Atom< LADA_TYPE >::frozen::X, "Could not reload freeze.\n");
  LADA_DOASSERT_TYPE;
  LADA_DOASSERT_SITE2;

  return 0;
}

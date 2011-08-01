#include "LaDaConfig.h"

#include<iostream>
#include<string>
#include<vector>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/erase.hpp>

#include <opt/debug.h>
#include <load_n_save/load/load.h>
#include <load_n_save/xml/parser.h>
#include <load_n_save/xml/printer.h>
#include <load_n_save/save/save.h>

#include "../structure.h"

#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE std::string 
#  define LADA_TYPE_INIT1 "Au"
#  define LADA_TYPE_INIT2 "Pd"
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE std::vector<std::string>
#  define LADA_TYPE_INIT1 "Au"
#  define LADA_TYPE_INIT2 "Au", "Pd"
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  namespace lns = LaDa::load_n_save;
  LADA_TEST_STRUCTURE< LADA_TYPE > structure;

  structure.set_cell(-0.5,0.5,0.5)
                    (0.5,-0.5,0.5)
                    (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, LADA_TYPE_INIT1);
  structure.add_atom(0.25,0.25,0.25, LADA_TYPE_INIT2);


  // print to xml
  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(structure)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  std::cout << xmlstring << "\n";
// LADA_DOASSERT(xmlstring == LADA_XML, "Error in output.");
//
// atom.pos = math::rVector3d(0,0,0);
// atom.freeze = Atom< LADA_TYPE >::frozen::NONE;
// LADA_CLEAR_TYPE;
// other = lns::xml::parse( sstr.str() );
// lns::load::Load loader;
// bool result = loader( *other, lns::ext(atom) );
// LADA_DOASSERT((atom.pos - math::rVector3d(0.5,-0.5,0.5)).squaredNorm() < 1e-12, "Could not reload pos.\n");
// LADA_DOASSERT(atom.freeze == Atom< LADA_TYPE >::frozen::X, "Could not reload freeze.\n");
// LADA_DOASSERT_TYPE
  return 0;
}

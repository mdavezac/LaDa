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

#if LADA_TEST_STRUCTURE == 0
#  define LADA_STRUCTURE StructureData
#  define LADA_CALL 
#  define LADA_ATOMS .atoms
#elif LADA_TEST_STRUCTURE == 1
#  define LADA_STRUCTURE Structure
#  define LADA_CALL   ()
#  define LADA_ATOMS 
#endif
#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE std::string 
#  define LADA_TYPE_INIT1 "Au"
#  define LADA_TYPE_INIT2 "Pd"
#  define LADA_XML "<Structure name=\"\" energy=\"0\" weight=\"1\" freeze=\"none\" scale=\"1\">  <Cell r0=\"-0.5 0.5 0.5\" r1=\"0.5 -0.5 0.5\" r2=\"0.5 0.5 -0.5\"/>  <Atom pos=\"0 0 0\" type=\"Au\" site=\"0\" freeze=\"none\"/>  <Atom pos=\"0.25 0.25 0.25\" type=\"Pd\" site=\"-1\" freeze=\"none\"/></Structure>"
#  define LADA_FIRST_SIZE true
#  define LADA_FIRST_TYPE structure LADA_ATOMS [0]->type == "Au"
#  define LADA_SECOND_SIZE true
#  define LADA_SECOND_TYPE structure LADA_ATOMS [1]->type == "Pd"
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE std::vector<std::string>
#  define LADA_TYPE_INIT1 "Au"
#  define LADA_TYPE_INIT2 "Au", "Pd"
#  define LADA_XML "<Structure name=\"\" energy=\"0\" weight=\"1\" freeze=\"none\" scale=\"1\">  <Cell r0=\"-0.5 0.5 0.5\" r1=\"0.5 -0.5 0.5\" r2=\"0.5 0.5 -0.5\"/>  <Atom pos=\"0 0 0\" type=\"Au\" site=\"0\" freeze=\"none\"/>  <Atom pos=\"0.25 0.25 0.25\" type=\"Au Pd\" site=\"-1\" freeze=\"none\"/></Structure>"
#  define LADA_FIRST_SIZE structure LADA_ATOMS [0]->type.size() == 1
#  define LADA_FIRST_TYPE structure LADA_ATOMS [0]->type[0] == "Au"
#  define LADA_SECOND_SIZE structure LADA_ATOMS [1]->type.size() == 2
#  define LADA_SECOND_TYPE     structure LADA_ATOMS [1]->type[0] == "Au" \
                           and structure LADA_ATOMS [1]->type[1] == "Pd"
#elif LADA_TEST_INCTYPE == 2
#  define LADA_TYPE std::set<std::string>
#  define LADA_TYPE_INIT1 "Au"
#  define LADA_TYPE_INIT2 "Au", "Pd"
#  define LADA_XML "<Structure name=\"\" energy=\"0\" weight=\"1\" freeze=\"none\" scale=\"1\">  <Cell r0=\"-0.5 0.5 0.5\" r1=\"0.5 -0.5 0.5\" r2=\"0.5 0.5 -0.5\"/>  <Atom pos=\"0 0 0\" type=\"Au\" site=\"0\" freeze=\"none\"/>  <Atom pos=\"0.25 0.25 0.25\" type=\"Au Pd\" site=\"-1\" freeze=\"none\"/></Structure>"
  template<class T> LADA_TYPE create_set(T const &_t)
    { LADA_TYPE result; result.insert(_t); return result; }
  template<class T> LADA_TYPE create_set(T const &_t, T const &_t1)
    { LADA_TYPE result; result.insert(_t); result.insert(_t1); return result; }
#  define LADA_FIRST_SIZE structure LADA_ATOMS [0]->type.size() == 1
#  define LADA_FIRST_TYPE structure LADA_ATOMS [0]->type == create_set("Au")
#  define LADA_SECOND_SIZE structure LADA_ATOMS [1]->type.size() == 2
#  define LADA_SECOND_TYPE structure LADA_ATOMS [1]->type == create_set("Au", "Pd")
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  namespace lns = LaDa::load_n_save;
  LADA_STRUCTURE< LADA_TYPE > structure;

  structure.set_cell(-0.5,0.5,0.5)
                    (0.5,-0.5,0.5)
                    (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, LADA_TYPE_INIT1)
                    (0.25,0.25,0.25, LADA_TYPE_INIT2);
  structure[0].site() = 0;


  // print to xml
  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(structure)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML, "Error in output.");

  structure.set_cell(0,0,0)
                    (0,0,0)
                    (0,0,0);
  structure.energy LADA_CALL= 666e0;
  structure.weight LADA_CALL= 666e0;
  structure.scale  LADA_CALL= 666e0;
  structure.freeze LADA_CALL= frozenstr::XX; 
  structure LADA_ATOMS .clear();
  other = lns::xml::parse( sstr.str() );
  lns::load::Load loader;
  bool result = loader( *other, lns::ext(structure) );
  for(size_t i(0); i < 3; ++i)
    for(size_t j(0); j < 3; ++j)
      LADA_DOASSERT( std::abs(structure.cell LADA_CALL (i,j) - (i==j ? -0.5:0.5)) < 1e-12,
                     "Could not reload cell.\n" )
  LADA_DOASSERT( std::abs(structure.energy LADA_CALL ) < 1e-12, "Could not reload energy.\n")
  LADA_DOASSERT( std::abs(structure.weight LADA_CALL -1) < 1e-12, "Could not reload weight.\n")
  LADA_DOASSERT( std::abs(structure.scale LADA_CALL -1) < 1e-12, "Could not reload scale.\n")
  LADA_DOASSERT( structure.freeze LADA_CALL == frozenstr::NONE, "Could not reload frozen.\n")
  LADA_DOASSERT( structure LADA_ATOMS .size() == 2, "Could not reload two atoms.\n")
  LADA_DOASSERT( (structure LADA_ATOMS [0].pos()).squaredNorm() < 1e-12,
                 "Could not reload first position.\n") 
  LADA_DOASSERT( (structure LADA_ATOMS [1].pos() - math::rVector3d(0.25,0.25,0.25)).squaredNorm() < 1e-12,
                 "Could not reload second position.\n") 
  LADA_DOASSERT( LADA_FIRST_SIZE, "Could not reload all first species.\n") 
  LADA_DOASSERT( LADA_FIRST_TYPE, "Could not reload first type.\n") 
  LADA_DOASSERT( LADA_SECOND_SIZE, "Could not reload all second species.\n") 
  LADA_DOASSERT( LADA_SECOND_TYPE, "Could not reload first type.\n") 
  LADA_DOASSERT( structure LADA_ATOMS[0].site() == 0, "Could not reload first site index.\n") 
  LADA_DOASSERT( structure LADA_ATOMS[1].site() == -1, "Could not reload second site index.\n") 
  return 0;
}

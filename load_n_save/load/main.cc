#include "LaDaConfig.h"

#define __PROGNAME__ "loadnsave"

#include "load.h"
#include "crystal/atom.h"
#include "crystal/structure.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../xml/printer.h"
#include "../action/enum.h"

#include "../save/save.h"

namespace lns = LaDa :: load_n_save;

int main(int argc, char *argv[]) 
{
  
  LaDa::Crystal::TStructure<std::string> structure;
  LaDa::Crystal::StrAtom atom;
  atom.pos[0] = 5;
  std::string string
    = "<Structure scale=\"5.45\">\n"
      "  <Cell a0=\"1 1 0\"\n"
      "        a1=\"0 1 0\"\n"
      "        a2=\"0 0 1\"/>\n"
      "  <Atom pos=\"0 0  3.1416\" type=Pd freeze=\"xyz\" />\n"
      "  <Atom pos=\"0 0  -3.1416\" type=Au freeze=\"xyz\" />\n"
      "</Structure>\n";
//   = "<Structure>\n"
//     "  <Cell row=\"1 0 0\" \n" 
//     "        row=\"0 1 0\"  \n" 
//     "        row=\"0 0 1\"/>\n" 
//     "</Structure>\n";

// boost::shared_ptr< LaDa::Crystal::Structure::t_Atom > a(new LaDa::Crystal::Structure::t_Atom);
// lns::ext(a);
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( string ) );
  std::cout << (bool) xml << "\n";
  if( not (bool) xml ) return 0;
      
  lns::load::Load loader;
  bool result = loader( *xml, lns::ext(structure) );
  std::cout << "result: " << result << "\n" << structure << "\n";
  

  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(structure)));
  std::cout << "\n ? \n";
  lns::xml::print(std::cout, *other);
  std::cout << "\n ? \n";


  return 0;
}

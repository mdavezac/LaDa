//
//  Version: $Id: main.cc 1231 2009-07-17 05:12:39Z davezac $
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define __PROGNAME__ "loadnsave"

#include <boost/mpl/at.hpp>

#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/adapted/mpl.hpp>

#include <boost/proto/proto.hpp>
#include <boost/proto/debug.hpp>
#include <boost/proto/proto_typeof.hpp>
#include <boost/proto/matches.hpp>
#include <boost/proto/make_expr.hpp>
#include <boost/proto/literal.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <math/serialize.h>
#include <crystal/atom.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <opt/bpo_macros.h>

#include "../grammar/grammar.h"
#include "../transforms/getsections.h"
#include "../transforms/getoptions.h"
#include "../transforms/getname.h"
#include "../transforms/gettag.h"
#include "../transforms/getaction.h"
#include "../transforms/getcontent.h"
#include "../transforms/gethelp.h"
#include "../xml/parser.h"
#include "../tree/tree.h"
#include "initializer.h"
#include "../push_back.h"

namespace LaDa
{
  namespace Crystal
  {
    template<class T_TYPE> template<class T_ARCHIVE, class TYPE_CONVERTER>
      bool TStructure<T_TYPE> :: lns_access_( T_ARCHIVE &_ar, TYPE_CONVERTER &_type )
      {
        namespace l = LaDa::load_n_save;
        namespace g = LaDa::load_n_save::grammar;
        namespace tgs = LaDa::load_n_save::tags;
        
        int const required(tgs::option::required);
        std::map<std::string, types::t_unsigned> freeze_x, freeze_y, freeze_z;
        freeze_x["x"] = FREEZE_XX;
        freeze_x["y"] = FREEZE_XY;
        freeze_x["z"] = FREEZE_XZ;
        freeze_y["all"] = FREEZE_XX | FREEZE_XY | FREEZE_XZ;
        freeze_y["x"] = FREEZE_XY;
        freeze_y["y"] = FREEZE_YY;
        freeze_y["z"] = FREEZE_YZ;
        freeze_y["all"] = FREEZE_XY | FREEZE_YY | FREEZE_YZ;
        freeze_z["x"] = FREEZE_XZ;
        freeze_z["y"] = FREEZE_YZ;
        freeze_z["z"] = FREEZE_ZZ;
        freeze_z["all"] = FREEZE_XZ | FREEZE_YZ | FREEZE_ZZ;
        if( _ar.is_loading() ) freeze = 0;
 
        bool loaded_cell =  _ar &
          ( 
            g::section( "Structure", g::tags=required )
              = (
                    g::option( "scale", g::action=scale, g::default_=2 )
                  + g::option( "name", g::action=name, g::default_=std::string() )
                  + g::option( "energy", g::action=energy, g::default_=666.666 )
                  + g::option( "weight", g::action=weight, g::default_=1.0 )
                  + (g::section("Cell", g::tags=required)
                      =(
                           (g::section("row", g::tags=required)
                             =(
                                  g::option( "x", g::tags=required, g::action=cell(0,0) )
                                + g::option( "y", g::tags=required, g::action=cell(0,1) )
                                + g::option( "z", g::tags=required, g::action=cell(0,2) )
                                + g::option
                                  (
                                    "freeze", 
                                    g::action=load_n_save::bitwise_map(freeze, freeze_x)
                                  )
                              ))
                         + (g::section("row", g::tags=required)
                             =(
                                  g::option( "x", g::tags=required, g::action=cell(1,0) )
                                + g::option( "y", g::tags=required, g::action=cell(1,1) )
                                + g::option( "z", g::tags=required, g::action=cell(1,2) )
                                + g::option
                                  (
                                    "freeze", 
                                    g::action=load_n_save::bitwise_map(freeze, freeze_y)
                                  )
                              ))
                         + (g::section("row", g::tags=required)
                             =(
                                  g::option( "x", g::tags=required, g::action=cell(2,0) )
                                + g::option( "y", g::tags=required, g::action=cell(2,1) )
                                + g::option( "z", g::tags=required, g::action=cell(2,2) )
                                + g::option
                                  (
                                    "freeze", 
                                    g::action=load_n_save::bitwise_map(freeze, freeze_z)
                                  )
                              ))
                       ))
                )  
          );
        if( not loaded_cell ) return false;
        t_Atom atom;
        bool loaded_atoms = _ar & ( g::section("Whatever") = g::ext(atom) ); 
        if( not loaded_atoms ) return false;
      }
  }
}





int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bf = boost::fusion;
  namespace proto = boost::proto;
  namespace lns = LaDa :: load_n_save;
  namespace lnsg = LaDa :: load_n_save :: grammar;
  using namespace lnsg;
  
  __TRYBEGIN
  __BPO_START__;
  __BPO_HIDDEN__;
  po::options_description all; 
  all.add(generic);
  po::options_description allnhidden;
  allnhidden.add(all).add(hidden);
  po::positional_options_description p; 
  p.add("input", 1); 
  po::variables_map vm; 
  po::store(po::command_line_parser(argc, argv). 
            options(allnhidden).positional(p).run(), vm); 
  po::notify(vm); 
  __BPO_PROGNAME__
  __BPO_VERSION__
  if( vm.count("help") )
  {
    std::cout << argv[0] << " is for testing only.\n" 
                 "Usage: " << argv[0] << " [options] file.xml\n" 
              << "  file.xml is an optional filename for XML input.\n" 
              << "  Default input is input.xml.\n\n" 
              << all << "\n";
    return 0; 
  }

  fs::path input( vm["input"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
  

  boost::shared_ptr<LaDa::Crystal::Lattice> lattice( LaDa::Crystal::read_lattice(input) );
  std::string string
    = "<Structure scale=\"5.45\" a = 5 >"                           "\n"
      "  <Cell>"                                              "\n"
      "    <row x=1 y=0 z =0 />"                              "\n"
      "    <row x=0 y=1 z =0 />"                              "\n"
      "    <row x=0 y=0 z =1 />"                              "\n"
      "  </Cell>"                                             "\n"
      "</Structure>"                                          "\n";
      
  if( not lattice ) return 1;
  LaDa::Crystal::Structure structure;
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( string ) );
  std::cout << (bool) xml << "\n";
  if( not (bool) xml ) return 0;
    
  lns::initializer::Initializer init;
  LaDa::Crystal::Atom atom;
  bool result = init( *xml, lnsg::ext(structure) );
  std::cout << "result: " << result << "\n";
  std::cout << structure << "\n";
  std::cout << structure.freeze << "\n";

  __BPO_CATCH__()

  return 0;
}

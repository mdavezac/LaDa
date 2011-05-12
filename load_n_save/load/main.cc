//
//  Version: $Id: main.cc 1293 2009-09-08 05:51:37Z davezac $
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define __PROGNAME__ "loadnsave"

#include "load.h"
#include "crystal/atom.h"
#include "../tags.h"
#include "../xpr/utilities.h"
#include "../xml/parser.h"
#include "../action/enum.h"

namespace lns = LaDa :: load_n_save;
struct A
{
  size_t a;
  std::string b;
  int c;
  template< class T_ARCHIVE >
    bool lns_access( T_ARCHIVE const &_ar )
    {
      std::map<std::string, LaDa::types::t_int> freeze_map;
      freeze_map["x"] = 1;
      freeze_map["y"] = 2;
      freeze_map["z"] = 4;
      freeze_map["t"] = 8;
      freeze_map["all"] = 15;

      return _ar &
        (
          lns::section("A") 
            << (
                     (
                          lns::option("a", lns::tag=lns::required, lns::action=a)
                       || lns::option("aa", lns::tag=lns::required, lns::action=a)
                     )
                  && lns::option("b", lns::action=b, lns::default_="dfsads")
                  && lns::option("c", lns::action=lns::enum_(c, freeze_map), lns::default_=0 )
               )
        );
    };
};

// namespace LaDa
// {
//   namespace load_n_save
//   {
//     template< class T_DATA >
//       struct lns_access
//       {
//         template< class T_ARCHIVE >
//           bool operator()( T_ARCHIVE const &_ar, T_DATA &_data ) const
//           {
//             return _ar &
//               (
//                 lns::section("A") 
//                   << (
//                           lns::option("a", lns::tag=lns::required, lns::action=_data.a)
//                         + lns::option("b", lns::action=_data.b)
//                         + lns::option("c", lns::action=_data.c)
//                      )
//               );
//           }
//       };
//   }
// }

int main(int argc, char *argv[]) 
{
  
  LaDa::Crystal::Atom atom;
  A a;
// lns::xpr::Section sec = lns::section( "Atom")
//                         << (
//                                lns::option( "x", lns::action=x )
//                              + lns::option( "y", lns::action=y )
//                              + lns::option( "z", lns::default_=1e-5, lns::action=z )
//                            );
  std::string string
    = "<A aa=4 b=\"3e-5\" c=t />\n";
//  = "<Atom x=4 y=0 z=\"-3.1416\" freeze=\"xyz\" />\n";

  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( string ) );
  std::cout << (bool) xml << "\n";
  if( not (bool) xml ) return 0;
      
  lns::load::Load loader;
  bool result = loader( *xml, lns::ext(a) );
// std::cout << "result: " << atom << "\n";
  std::cout << "result: " << result << " " << a.a << " " << a.b << " " << a.c << "\n";


  return 0;
}

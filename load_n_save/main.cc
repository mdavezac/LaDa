//
//  Version: $Id: main.cc 1200 2009-06-22 05:18:21Z davezac $
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

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


#include "grammar/grammar.h"
#include "transforms/getsections.h"
#include "transforms/getoptions.h"
#include "transforms/getname.h"
#include "transforms/getaction.h"
#include "transforms/getcontent.h"
#include "initializer/initializer.h"
#include <crystal/atom.h>

namespace LaDa
{
  namespace load_n_save
  {
    template<class T_ARCHIVE>
      bool lns_access( T_ARCHIVE const& _ar, LaDa::Crystal::Atom &_atom )
      {
        return _ar &
          ( 
            grammar::section( "Atom" )
              = (
                    grammar::option("x")[_atom.x]
                  + grammar::option("y")[_atom.y]
                  + grammar::option("z")[_atom.z]
                )
          );
      }
  }
}

int main(int argc, char *argv[]) 
{
  namespace bf = boost::fusion;
  namespace proto = boost::proto;
  namespace lns = LaDa :: load_n_save;
  namespace lnsg = LaDa :: load_n_save :: grammar;
      
  std::string string = "here";
  LaDa::types::t_unsigned uint = 5u;
  LaDa::types::t_int integer = 6;
  LaDa::types::t_real real = 6.3;
# define tree \
   (\
     lnsg::section("I") \
       = ( \
           lnsg::option("a") = ( lnsg::lit("whatever") == lnsg::lit(5) ) || lnsg::lit(6e0) \
         )[ lnsg::var( string ) ]\
   )\
   + ( \
       lnsg::section("II") = \
       ( \
           ( ( lnsg::option("A") = 5e0 ) ) + lnsg::option("B")\
         + lnsg::section("sub i")[lnsg::call(integer) ] \
         + ( lnsg::section("sub ii") = (   lnsg::option("b")[lnsg::var(uint)]\
                                         + lnsg::option("c") ) ) \
       )\
     )\
   + lnsg::section("III")[ lnsg::call( integer ) ]
// proto::display_expr
// (  
//   tree
// );
  boost::fusion::for_each
  ( 
    lns::transform::GetSections()( tree ), 
    print_section() 
  );
//
// lns::Print()( lns::transform::Sections()( tree, boost::fusion::nil() ) );

    

  return 0;
}

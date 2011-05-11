//
//  Version: $Id: main.cc 1167 2009-06-07 23:22:42Z davezac $
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>

#include "parser.h"
#include "printer.h"

int main(int argc, char *argv[]) 
{
  std::string string = "  \n<Whatever option=\"a\"> <Text u= \"5ag\" \n b=36 d=5 /> \n</Whatever>";

  namespace lns = LaDa::load_n_save;
  boost::shared_ptr< lns::tree::Base > xml( lns::xml::parse( string ) );
  std::cout << (bool) xml << "\n";
  if( not (bool) xml ) return 0;
  lns::xml::print( std::cout, *xml );
    

  return 0;
}

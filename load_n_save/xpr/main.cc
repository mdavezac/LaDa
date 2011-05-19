#include "LaDaConfig.h"

#include "utilities.h"
#include "section.h"


int main(int argc, char *argv[]) 
{
  namespace x = LaDa::load_n_save;
  x::xpr::Section b = ( x::section("First")  );
  x::xpr::Section c = (
                                 x::option( "Whatever", x::help="not really" )
                              || x::option( "Other", x::help="not really", x::tag=5  ) 
                              && x::option( "N", x::help="not really", x::tag=5  )
                              && x::option( "c" ) 
                              && x::option( "d" )
                      );
  x::xpr::Section a = (x::section("First")) << c;
  std::cout << a.sequence() << "\n";
  std::cout << a.print() << "\n";
  return 0;
}

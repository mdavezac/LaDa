#include "motu.h"
using LaDa::MotU;

int main(int argc, char *argv[]) 
{
  try
  {
    MotU motu;
    if ( not  motu.Load("input.xml") )
      throw ""; 
    motu.run();
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running lada" << std::endl
              << e.what();
  }
  return 0;
}
